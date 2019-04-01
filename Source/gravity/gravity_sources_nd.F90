module gravity_sources_module

  implicit none

  public

contains

  subroutine ca_gsrc(lo, hi, &
                     domlo, domhi, &
                     uold, uold_lo, uold_hi, &
#ifdef SELF_GRAVITY
                     phi, phi_lo, phi_hi, &
                     grav, grav_lo, grav_hi, &
#endif
                     source, src_lo, src_hi, &
                     dx, dt, time) bind(C, name="ca_gsrc")

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO, HALF, ONE
#ifndef AMREX_USE_CUDA
    use amrex_error_module, only: amrex_error
#endif
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN
    use castro_util_module, only: position ! function
    use prob_params_module, only: center
#ifdef HYBRID_MOMENTUM
    use meth_params_module, only: UMR, UMP
    use hybrid_advection_module, only: set_hybrid_momentum_source
#endif
#ifndef SELF_GRAVITY
    use meth_params_module, only: const_grav
    use prob_params_module, only: dim
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uold_lo(3), uold_hi(3)
#ifdef SELF_GRAVITY
    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
#endif
    integer,  intent(in   ) :: src_lo(3), src_hi(3)

    real(rt), intent(in   ) :: uold(uold_lo(1):uold_hi(1),uold_lo(2):uold_hi(2),uold_lo(3):uold_hi(3),NVAR)
#ifdef SELF_GRAVITY
    real(rt), intent(in   ) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(rt), intent(in   ) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3),3)
#endif
    real(rt), intent(inout) :: source(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt, time

    real(rt) :: rho, rhoInv
    real(rt) :: Sr(3), SrE
    real(rt) :: old_ke, new_ke
    real(rt) :: loc(3)

    integer  :: i, j, k

    ! Temporary array for holding the update to the state.
    
    real(rt) :: src(NVAR)

    ! Temporary array for seeing what the new state would be if the update were applied here.

    real(rt) :: snew(NVAR)

    !$gpu

    ! Initialize the update and temporary state to zero. We only need to do this once outside
    ! the loop, since the array access pattern is consistent across loop iterations.

    Sr(:) = ZERO
    src(:) = ZERO
    snew(:) = ZERO

    ! For constant gravity, we can just initialize the gravitational acceleration here.

#ifndef SELF_GRAVITY
    Sr(dim) = const_grav
#endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rho    = uold(i,j,k,URHO)
             rhoInv = ONE / rho

             loc = position(i,j,k) - center

             src = ZERO
             snew = uold(i,j,k,:)

             old_ke = HALF * sum(snew(UMX:UMZ)**2) * rhoInv

#ifdef SELF_GRAVITY
             Sr = rho * grav(i,j,k,:)
#else
             Sr(dim) = rho * const_grav
#endif

             src(UMX:UMZ) = Sr

             snew(UMX:UMZ) = snew(UMX:UMZ) + dt * src(UMX:UMZ)

#ifdef HYBRID_MOMENTUM
             call set_hybrid_momentum_source(loc, src(UMR:UMP), Sr)

             snew(UMR:UMP) = snew(UMR:UMP) + dt * src(UMR:UMP)
#endif

             ! The conservative energy formulation does not strictly require
             ! any energy source-term here, because it depends only on the
             ! fluid motions from the hydrodynamical fluxes which we will only
             ! have when we get to the 'corrector' step. Nevertheless we add a
             ! predictor energy source term in the way that the other methods
             ! do, for consistency. We will fully subtract this predictor value
             ! during the corrector step, so that the final result is correct.

             SrE = dot_product(uold(i,j,k,UMX:UMZ) * rhoInv, Sr)

             src(UEDEN) = SrE

             snew(UEDEN) = snew(UEDEN) + dt * SrE

             ! Add to the outgoing source array.

             source(i,j,k,:) = source(i,j,k,:) + src

          enddo
       enddo
    enddo

  end subroutine ca_gsrc

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::

  subroutine ca_corrgsrc(lo, hi, &
                         domlo, domhi, &
                         uold, uo_lo, uo_hi, &
                         unew, un_lo, un_hi, &
#ifdef SELF_GRAVITY
                         gold, go_lo, go_hi, &
                         gnew, gn_lo, gn_hi, &
                         gxold, xo_lo, xo_hi, &
                         gxnew, xn_lo, xn_hi, &
                         gyold, yo_lo, yo_hi, &
                         gynew, yn_lo, yn_hi, &
                         gzold, zo_lo, zo_hi, &
                         gznew, zn_lo, zn_hi, &
#endif
                         vol, vol_lo, vol_hi, &
                         flux1, f1_lo, f1_hi, &
                         flux2, f2_lo, f2_hi, &
                         flux3, f3_lo, f3_hi, &
                         source, sr_lo, sr_hi, &
                         dx, dt, time) bind(C, name="ca_corrgsrc")

    use amrex_fort_module, only: rt => amrex_real
#ifndef AMREX_USE_CUDA
    use amrex_error_module, only: amrex_error
#endif
    use amrex_constants_module, only: ZERO, HALF, ONE, TWO
    use amrex_mempool_module, only: bl_allocate, bl_deallocate
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, &
                                  gravity_type_int, get_g_from_phi
    use prob_params_module, only: dg, center, physbc_lo, physbc_hi, Symmetry
    use castro_util_module, only: position ! function
#ifdef HYBRID_MOMENTUM
    use meth_params_module, only: UMR, UMP
    use hybrid_advection_module, only: set_hybrid_momentum_source
#endif
#ifndef SELF_GRAVITY
    use meth_params_module, only: const_grav
    use prob_params_module, only: dim
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uo_lo(3), uo_hi(3)
    integer,  intent(in   ) :: un_lo(3), un_hi(3)
#ifdef SELF_GRAVITY
    integer,  intent(in   ) :: go_lo(3), go_hi(3)
    integer,  intent(in   ) :: gn_lo(3), gn_hi(3)
    integer,  intent(in   ) :: xo_lo(3), xo_hi(3)
    integer,  intent(in   ) :: xn_lo(3), xn_hi(3)
    integer,  intent(in   ) :: yo_lo(3), yo_hi(3)
    integer,  intent(in   ) :: yn_lo(3), yn_hi(3)
    integer,  intent(in   ) :: zo_lo(3), zo_hi(3)
    integer,  intent(in   ) :: zn_lo(3), zn_hi(3)
#endif
    integer,  intent(in   ) :: vol_lo(3), vol_hi(3)
    integer,  intent(in   ) :: f1_lo(3), f1_hi(3)
    integer,  intent(in   ) :: f2_lo(3), f2_hi(3)
    integer,  intent(in   ) :: f3_lo(3), f3_hi(3)

    integer,  intent(in   ) :: sr_lo(3), sr_hi(3)

    ! Old and new time state data

    real(rt), intent(in   ) :: uold(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),NVAR)
    real(rt), intent(in   ) :: unew(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3),NVAR)

#ifdef SELF_GRAVITY
    ! Old and new time gravitational acceleration, cell-centered

    real(rt), intent(in   ) :: gold(go_lo(1):go_hi(1),go_lo(2):go_hi(2),go_lo(3):go_hi(3),3)
    real(rt), intent(in   ) :: gnew(gn_lo(1):gn_hi(1),gn_lo(2):gn_hi(2),gn_lo(3):gn_hi(3),3)

    ! Old and new time gravitational acceleration, edge-centered

    real(rt), intent(in   ) :: gxold(xo_lo(1):xo_hi(1),xo_lo(2):xo_hi(2),xo_lo(3):xo_hi(3))
    real(rt), intent(in   ) :: gxnew(xn_lo(1):xn_hi(1),xn_lo(2):xn_hi(2),xn_lo(3):xn_hi(3))
    real(rt), intent(in   ) :: gyold(yo_lo(1):yo_hi(1),yo_lo(2):yo_hi(2),yo_lo(3):yo_hi(3))
    real(rt), intent(in   ) :: gynew(yn_lo(1):yn_hi(1),yn_lo(2):yn_hi(2),yn_lo(3):yn_hi(3))
    real(rt), intent(in   ) :: gzold(zo_lo(1):zo_hi(1),zo_lo(2):zo_hi(2),zo_lo(3):zo_hi(3))
    real(rt), intent(in   ) :: gznew(zn_lo(1):zn_hi(1),zn_lo(2):zn_hi(2),zn_lo(3):zn_hi(3))
#endif

    ! Cell volume

    real(rt), intent(in   ) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))

    ! Hydrodynamical mass fluxes

    real(rt), intent(in   ) :: flux1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),f1_lo(3):f1_hi(3))
    real(rt), intent(in   ) :: flux2(f2_lo(1):f2_hi(1),f2_lo(2):f2_hi(2),f2_lo(3):f2_hi(3))
    real(rt), intent(in   ) :: flux3(f3_lo(1):f3_hi(1),f3_lo(2):f3_hi(2),f3_lo(3):f3_hi(3))

    ! The source term to send back

    real(rt), intent(inout) :: source(sr_lo(1):sr_hi(1),sr_lo(2):sr_hi(2),sr_lo(3):sr_hi(3),NVAR)

    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt, time

    integer  :: i, j, k

    real(rt) :: Sr_old(3), Sr_new(3), Srcorr(3)
    real(rt) :: vold(3), vnew(3)
    real(rt) :: SrE_old, SrE_new, SrEcorr
    real(rt) :: rhoo, rhooinv, rhon, rhoninv

    real(rt) :: old_ke, new_ke
    real(rt) :: loc(3)

    real(rt) :: hdtInv

    real(rt) :: phi, phixl, phixr, phiyl, phiyr, phizl, phizr
    real(rt) :: g(3), gl(3), gr(3)

    real(rt) :: src(NVAR)

    ! Temporary array for seeing what the new state would be if the update were applied here.

    real(rt) :: snew(NVAR)

    !$gpu

    Sr_old(:) = ZERO
    Sr_new(:) = ZERO
    Srcorr(:) = ZERO
    src(:) = ZERO
    snew(:) = ZERO

    hdtInv = HALF / dt

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             rhoo    = uold(i,j,k,URHO)
             rhooinv = ONE / uold(i,j,k,URHO)

             rhon    = unew(i,j,k,URHO)
             rhoninv = ONE / unew(i,j,k,URHO)

             snew = unew(i,j,k,:)

             old_ke = HALF * sum(snew(UMX:UMZ)**2) * rhoninv

             ! Define old source terms

             vold = uold(i,j,k,UMX:UMZ) * rhooinv

#ifdef SELF_GRAVITY
             Sr_old = rhoo * gold(i,j,k,:)
#else
             Sr_old(dim) = rhoo * const_grav
#endif
             SrE_old = dot_product(vold, Sr_old)

             ! Define new source terms

             vnew = snew(UMX:UMZ) * rhoninv

#ifdef SELF_GRAVITY
             Sr_new = rhon * gnew(i,j,k,:)
#else
             Sr_new(dim) = rhon * const_grav
#endif
             SrE_new = dot_product(vnew, Sr_new)

             ! Define corrections to source terms

             Srcorr = HALF * (Sr_new - Sr_old)

             ! Correct momenta

             src(UMX:UMZ) = Srcorr

             snew(UMX:UMZ) = snew(UMX:UMZ) + dt * src(UMX:UMZ)

#ifdef HYBRID_MOMENTUM
             call set_hybrid_momentum_source(loc, src(UMR:UMP), Srcorr)

             snew(UMR:UMP) = snew(UMR:UMP) + dt * src(UMR:UMP)
#endif

             ! Correct energy

             ! First, subtract the predictor step we applied earlier.

             SrEcorr = - SrE_old

             ! For an explanation of this approach, see wdmerger paper I.
             ! The main idea is that we are evaluating the change of the
             ! potential energy at zone edges and applying that in an equal
             ! and opposite sense to the gas energy. The physics is described
             ! in Section 2.4. We use a slightly different formulation than
             ! listed in the paper, which relies on the (second-order) property:
             ! g_{i+1/2} = -( phi_{i+1} - phi_{i} ) / dx.

#ifdef SELF_GRAVITY

             ! Construct the time-averaged edge-centered gravity.
             ! Note that what comes in grad(phi) so we need to negate it.

             gl(1) = -HALF * (gxold(i,        j,k) + gxnew(i,        j,k))
             gr(1) = -HALF * (gxold(i+1*dg(1),j,k) + gxnew(i+1*dg(1),j,k))

             gl(2) = -HALF * (gyold(i,j        ,k) + gynew(i,j,        k))
             gr(2) = -HALF * (gyold(i,j+1*dg(2),k) + gynew(i,j+1*dg(2),k))

             gl(3) = -HALF * (gzold(i,j,k        ) + gznew(i,j,k        ))
             gr(3) = -HALF * (gzold(i,j,k+1*dg(3)) + gznew(i,j,k+1*dg(3)))

#else

             gl(:) = ZERO
             gr(:) = ZERO

             gl(dim) = const_grav
             gr(dim) = const_grav

#endif

             SrEcorr = SrEcorr + hdtInv * ( flux1(i        ,j,k) * gl(1) * dx(1) + &
                                            flux1(i+1*dg(1),j,k) * gr(1) * dx(1) + &
                                            flux2(i,j        ,k) * gl(2) * dx(2) + &
                                            flux2(i,j+1*dg(2),k) * gr(2) * dx(2) + &
                                            flux3(i,j,k        ) * gl(3) * dx(3) + &
                                            flux3(i,j,k+1*dg(3)) * gr(3) * dx(3) ) / vol(i,j,k)

             src(UEDEN) = SrEcorr

             snew(UEDEN) = snew(UEDEN) + dt * SrEcorr

             ! Add to the outgoing source array.

             source(i,j,k,:) = source(i,j,k,:) + src

          enddo
       enddo
    enddo

  end subroutine ca_corrgsrc

end module gravity_sources_module
