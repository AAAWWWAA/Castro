#ifndef _Castro_advection_util_H_
#define _Castro_advection_util_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_HOST_DEVICE
inline
void normalize_species_fluxes(int i, int j, int k, amrex::Array4<amrex::Real> flux)
{

    // Normalize the fluxes of the mass fractions so that
    // they sum to 0.  This is essentially the CMA procedure that is
    // defined in Plewa & Muller, 1999, A&A, 342, 179.

    amrex::Real sum = 0.0;
    amrex::Real fac = 1.0;

    for (int n = 7; n < 8; ++n) {
        sum = sum + flux(i,j,k,n);
    }

    if (sum != 0.0) {
        fac = flux(i,j,k,0) / sum;
    }

    for (int n = 7; n < 8; ++n) {
        flux(i,j,k,n) *= fac;
    }

}

#endif
