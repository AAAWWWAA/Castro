#include "Castro.H"
#include "Castro_F.H"

#include "Gravity.H"

using namespace amrex;

#ifdef GRAVITY
#ifdef ROTATION
void Castro::scf_relaxation() {

    int finest_level = parent->finestLevel();

    MultiFab& S_new = get_new_data(State_Type);

    int j = 1;
    int relax_max_iterations = 30;

    const Real* dx = parent->Geom(level).CellSize();

    // First do some sanity checks.

    // Maximum density must be set.

    if (scf_maximum_density <= 0) {
        amrex::Error("castro.scf_maximum_density must be set for SCF relaxation.");
    }

    // Equatorial radius and polar radius must both
    // be positive, and the polar radius cannot be
    // larger than the equatorial radius.

    if (scf_equatorial_radius <= 0) {
        amrex::Error("Equatorial radius must be positive for SCF relaxation.");
    }

    if (scf_polar_radius <= 0) {
        amrex::Error("Polar radius must be positive for SCF relaxation.");
    }

    if (scf_polar_radius > scf_equatorial_radius) {
        amrex::Error("Polar radius must not be larger than equatorial radius for SCF relaxation.");
    }

    // Grab the value for the solar mass, we'll need it later.

    Real M_solar;
    scf_get_solar_mass(&M_solar);

    Real time = state[State_Type].curTime();

    // Do the initial relaxation setup.

    scf_setup_relaxation();

    // To do the relaxation loop, we need to know the target maximum
    // enthalpy on the domain. So we'll loop through the grid and
    // calculate what it should be, based on the density, temperature,
    // and composition. We do not specify the latter two -- we assume
    // that these are filled as desired in initdata. So, to get the
    // maximum enthalpy, we'll loop through the grid, using the composition
    // and temperature in each zone, combined with the requested maximum
    // density, call the EOS to calculate the enthalpy for that combination,
    // and then take the maximum of that. (Note that the current density on
    // the grid is treated as an initial guess and does not need to reach
    // the requested maximum.) The SCF procedure is only intended for a setup
    // with uniform temperature and composition, so in that case this procedure
    // should yield the same result as just picking a zone at random and using
    // its temperature and composition. The reason for doing it this way is so
    // that the user doesn't have to request temperature and composition in their
    // inputs file (which is messy in general) -- they can tell us what composition
    // they want by filling it in the initialization.

    Real target_h_max = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(max:target_h_max)
#endif
    for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        scf_calculate_target_h_max(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                                   BL_TO_FORTRAN_ANYD(S_new[mfi]),
                                   &target_h_max);

    }

    ParallelDescriptor::ReduceRealMax(target_h_max);

    // Get the phi and phi_rot MultiFabs.

    MultiFab& phi = get_new_data(PhiGrav_Type);
    MultiFab& phi_rot = get_new_data(Rotation_Type);

    // Construct a local MultiFab for the rotational psi.
    // This will not change over the loop iterations.

    MultiFab psi(grids, dmap, 1, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(psi,true); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        ca_fill_rotational_psi(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                               BL_TO_FORTRAN_ANYD(psi[mfi]),
                               AMREX_ZFILL(dx), time);

    }

    // Iterate until the system is relaxed by filling the level data
    // and then doing a multilevel gravity solve.

    int is_relaxed = 0;

    while (j <= relax_max_iterations) {

        // First step is to find the rotational frequency.

        Real phi_A = 0.0;
        Real psi_A = 0.0;
        Real phi_B = 0.0;
        Real psi_B = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:phi_A, psi_A, phi_B, psi_B)
#endif
        for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

            scf_update_for_omegasq(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                                   BL_TO_FORTRAN_ANYD(S_new[mfi]),
                                   BL_TO_FORTRAN_ANYD(phi[mfi]),
                                   BL_TO_FORTRAN_ANYD(psi[mfi]),
                                   AMREX_ZFILL(dx),
                                   &phi_A, &psi_A, &phi_B, &psi_B);

        }

        ParallelDescriptor::ReduceRealSum(phi_A);
        ParallelDescriptor::ReduceRealSum(psi_A);
        ParallelDescriptor::ReduceRealSum(phi_B);
        ParallelDescriptor::ReduceRealSum(psi_B);

        // Now update the square of the rotation frequency, following Hachisu (Equation 16).
        // Deal carefully with the special case where phi_A and phi_B are equal -- we assume
        // this is the non-rotating case, and we want to handle that gracefully. Since the
        // rotation module knows to treat rotational_period = 0 as effectively disabling
        // rotation, we'll use that approach.

        if (std::abs(phi_A - phi_B) / std::abs(phi_A) < 1.e-6) {

            rotational_period = 0.0;

        } else {

            Real omegasq = -(phi_A - phi_B) / (psi_A - psi_B);

            if (omegasq < 0.0 && ParallelDescriptor::IOProcessor()) {
                std::cerr << "Omega squared = " << omegasq << " is negative in the relaxation step; aborting." << std::endl;
                amrex::Error();
            }

            Real omega = sqrt(omegasq);

            // Rotational period is 2 pi / omega.

            rotational_period = 2.0 * M_PI / omega;

        }

        // Now save the updated rotational frequency in the Fortran module.

        set_rot_period(&rotational_period);

        // With the updated period, we can construct the updated rotational
        // potential, which will be used in the remaining steps below.

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(phi_rot, true); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

            ca_fill_rotational_potential(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                                         BL_TO_FORTRAN_ANYD(phi_rot[mfi]),
                                         dx, time);

        }



        // Second step is to evaluate the Bernoulli constant.

        Real bernoulli = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:bernoulli)
#endif
        for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

            scf_get_bernoulli_const(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                                    BL_TO_FORTRAN_ANYD(S_new[mfi]),
                                    BL_TO_FORTRAN_ANYD(phi[mfi]),
                                    BL_TO_FORTRAN_ANYD(phi_rot[mfi]),
                                    AMREX_ZFILL(dx), &bernoulli);

        }

        ParallelDescriptor::ReduceRealSum(bernoulli);



        // Third step is to construct the enthalpy field and
        // find the maximum enthalpy for the star.

        MultiFab enthalpy(grids, dmap, 1, 0);
        enthalpy.setVal(0.0);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

            scf_construct_enthalpy(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                                   BL_TO_FORTRAN_ANYD(S_new[mfi]),
                                   BL_TO_FORTRAN_ANYD(phi[mfi]),
                                   BL_TO_FORTRAN_ANYD(phi_rot[mfi]),
                                   BL_TO_FORTRAN_ANYD(enthalpy[mfi]),
                                   AMREX_ZFILL(dx), bernoulli);

        }

        Real actual_h_max = enthalpy.max(0);

        Real Linf_norm = 0.0;

        // Finally, update the density using the enthalpy field.

#ifdef _OPENMP
#pragma omp parallel reduction(max:Linf_norm)
#endif
        for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

            scf_update_density(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                               BL_TO_FORTRAN_ANYD(S_new[mfi]),
                               BL_TO_FORTRAN_ANYD(enthalpy[mfi]),
                               AMREX_ZFILL(dx), actual_h_max, target_h_max,
                               &Linf_norm);

        }

        ParallelDescriptor::ReduceRealMax(Linf_norm);

        // Update diagnostic quantities.

        Real kin_eng = 0.0;
        Real pot_eng = 0.0;
        Real int_eng = 0.0;
        Real mass = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:kin_eng,pot_eng,int_eng,mass)
#endif
        for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

            scf_diagnostics(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                            BL_TO_FORTRAN_ANYD(S_new[mfi]),
                            BL_TO_FORTRAN_ANYD(phi[mfi]),
                            BL_TO_FORTRAN_ANYD(phi_rot[mfi]),
                            AMREX_ZFILL(dx),
                            &kin_eng, &pot_eng, &int_eng, &mass);

        }

        ParallelDescriptor::ReduceRealSum(kin_eng);
        ParallelDescriptor::ReduceRealSum(pot_eng);
        ParallelDescriptor::ReduceRealSum(int_eng);
        ParallelDescriptor::ReduceRealSum(mass);

        Real virial_error = std::abs(2.0 * kin_eng + pot_eng + 3.0 * int_eng) / std::abs(pot_eng);

        // Now check to see if we're converged.

        int is_relaxed = 0;

        if (Linf_norm < scf_relax_tol) {
            is_relaxed = 1;
        }

        if (ParallelDescriptor::IOProcessor()) {

            std::cout << std::endl << std::endl;
            std::cout << "   Relaxation iterations completed: " << j << std::endl;
            std::cout << "   L-infinity norm of residual (relative to old state): " << Linf_norm << std::endl;
            std::cout << "   Rotational period (s): " << rotational_period << std::endl;
            std::cout << "   Kinetic energy: " << kin_eng << std::endl;
            std::cout << "   Potential energy: " << pot_eng << std::endl;
            std::cout << "   Internal energy: " << int_eng << std::endl;
            std::cout << "   Virial error: " << virial_error << std::endl;
            std::cout << "   Mass: " << mass / M_solar << " solar masses" << std::endl;

            if (is_relaxed == 1) {
                std::cout << "  Relaxation completed!" << std::endl;
            }

            std::cout << std::endl << std::endl;

        }

        //	for (int k = finest_level-1; k >= 0; k--)
        //	  getLevel(k).avgDown();

        gravity->multilevel_solve_for_new_phi(level,finest_level);

        if (is_relaxed == 1) break;

        j++;

    }

    for (int k = 0; k <= parent->finestLevel(); k++)
    {
        MultiFab& grav_new = getLevel(k).get_new_data(Gravity_Type);
        gravity->get_new_grav_vector(k,grav_new,time);
    }

}
#endif
#endif
