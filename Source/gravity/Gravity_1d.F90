module gravity_1D_module

  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_test_residual(lo, hi, &
       rhs, rhlo, rhhi,  &
       ecx, ecxlo, ecxhi, &
       dx, problo, coord_type) bind(C, name="ca_test_residual")

    use amrex_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer , intent(in   ) :: lo(3),hi(3)
    integer , intent(in   ) :: rhlo(3), rhhi(3)
    integer , intent(in   ) :: ecxlo(3), ecxhi(3)
    integer , intent(in   ) :: coord_type
    real(rt), intent(inout) :: rhs(rhlo(1):rhhi(1),rhlo(2):rhhi(2),rhlo(3):rhhi(3))
    real(rt), intent(in   ) :: ecx(ecxlo(1):ecxhi(1),ecxlo(2):ecxhi(2),ecxlo(3):ecxhi(3))
    real(rt), intent(in   ) :: dx(3),problo(3)

    real(rt)         :: lapphi
    real(rt)         :: rlo,rhi,rcen
    integer          :: i, j, k

    j = lo(2)
    k = lo(3)

    ! Cartesian
    if (coord_type .eq. 0) then

       do i=lo(1),hi(1)
          lapphi = (ecx(i+1,j,k)-ecx(i,j,k)) / dx(1)
          rhs(i,j,k) = rhs(i,j,k) - lapphi
       enddo

       ! r-z
    else if (coord_type .eq. 1) then

       do i=lo(1),hi(1)
          rlo  = problo(1) + dble(i)*dx(1)
          rhi  = rlo + dx(1)
          rcen = HALF * (rlo+rhi)
          lapphi = (rhi*ecx(i+1,j,k)-rlo*ecx(i,j,k)) / (rcen*dx(1))
          rhs(i,j,k) = rhs(i,j,k) - lapphi
       enddo

       ! spherical
    else if (coord_type .eq. 2) then

       do i=lo(1),hi(1)
          rlo  = problo(1) + dble(i)*dx(1)
          rhi  = rlo + dx(1)
          rcen = HALF * (rlo+rhi)
          lapphi = (rhi**2 * ecx(i+1,j,k)-rlo**2 * ecx(i,j,k)) / (rcen**2 * dx(1))
          rhs(i,j,k) = rhs(i,j,k) - lapphi
       enddo

    else
       print *,'Bogus coord_type in test_residual ' ,coord_type
       call amrex_error("Error:: Gravity_1d.f90 :: ca_test_residual")
    end if

  end subroutine ca_test_residual



  subroutine ca_compute_radial_mass (lo,hi,dx,dr,&
                                     state,r_lo,r_hi, &
                                     radial_mass,radial_vol,problo, &
                                     n1d,drdxfac,level) bind(C, name="ca_compute_radial_mass")

    use amrex_constants_module, only: HALF, FOUR3RD, M_PI
    use prob_params_module, only: center, Symmetry, physbc_lo, coord_type
    use meth_params_module, only: NVAR, URHO

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer , intent(in   ) :: lo(3),hi(3)
    integer , intent(in   ) :: r_lo(3),r_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), value, intent(in   ) :: dr
    real(rt), intent(in   ) :: problo(3)

    integer , value, intent(in   ) :: n1d,drdxfac,level
    real(rt), intent(inout) :: radial_mass(0:n1d-1)
    real(rt), intent(inout) :: radial_vol (0:n1d-1)

    real(rt), intent(in   ) :: state(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),NVAR)

    integer          :: i, j, k, index
    integer          :: ii
    real(rt)         :: r
    real(rt)         :: dx_frac, fac, vol
    real(rt)         :: lo_i, rlo, rhi

#ifndef AMREX_USE_CUDA
    if (physbc_lo(1) .ne. Symmetry) then
       call amrex_error("Error: Gravity_1d.f90 :: 1D gravity assumes symmetric lower boundary.")
    endif

    if (coord_type .ne. 2) then
       call amrex_error("Error: Gravity_1d.f90 :: 1D gravity assumes spherical coordinates.")
    endif
#endif

    fac = dble(drdxfac)
    dx_frac = dx(1) / fac
    j = lo(2)
    k = lo(3)

    do i = lo(1), hi(1)

       r = abs(problo(1) + (dble(i) + HALF) * dx(1) - center(1))

       index = int(r / dr)

       if (index .gt. n1d-1) then
#ifndef AMREX_USE_CUDA
          if (level .eq. 0) then
             print *,'   '
             print *,'>>> Error: Gravity_1d::ca_compute_radial_mass ', i
             print *,'>>> ... index too big: ', index,' > ',n1d-1
             print *,'>>> ... at i     : ', i
             print *,'    '
             call amrex_error("Error:: Gravity_1d.f90 :: ca_compute_radial_mass")
          end if
#endif
       else

          ! Note that we assume we are in spherical coordinates in 1D or we wouldn't be
          ! doing monopole gravity.

          lo_i = problo(1) + dble(i) * dx(1) - center(1)

          do ii = 0, drdxfac-1

             r   = abs(lo_i + (dble(ii  ) + HALF) * dx_frac)
             rlo = abs(lo_i +  dble(ii  )         * dx_frac)
             rhi = abs(lo_i +  dble(ii+1)         * dx_frac)

             vol = FOUR3RD * M_PI * (rhi**3 - rlo**3)

             index = int(r / dr)

             if (index .le. n1d-1) then
                radial_mass(index) = radial_mass(index) + vol * state(i,j,k,URHO)
                radial_vol (index) = radial_vol (index) + vol
             end if

          enddo

       endif

    enddo

  end subroutine ca_compute_radial_mass



  subroutine ca_put_radial_grav(lo,hi,dx,dr,&
                                grav,g_lo,g_hi, &
                                radial_grav,problo,n1d,level) bind(C, name="ca_put_radial_grav")

    use prob_params_module, only: center
    use amrex_constants_module, only: ZERO, HALF, TWO

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer , intent(in   ) :: lo(3),hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), value, intent(in   ) :: dr
    real(rt), intent(in   ) :: problo(3)

    integer , value, intent(in   ) :: n1d,level
    real(rt), intent(in   ) :: radial_grav(0:n1d-1)

    integer , intent(in   ) :: g_lo(3),g_hi(3)
    real(rt), intent(inout) :: grav(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),3)

    integer          :: i, j, k, index
    real(rt)         :: r, mag_grav
    real(rt)         :: cen, xi, slope, glo, gmd, ghi, minvar, maxvar

    !$gpu

    ! Note that we are interpolating onto the entire range of grav,
    ! including the ghost cells. Taking the absolute value of r ensures
    ! that we will get the correct behavior even for the ghost zones with
    ! negative indices, which have a reflecting boundary condition.

    j = lo(2)
    k = lo(3)

    grav(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:3) = ZERO

    do i = lo(1), hi(1)

       r = abs(problo(1) + (dble(i) + HALF) * dx(1) - center(1))

       index = int(r / dr)

       cen = (dble(index) + HALF) * dr
       xi = r - cen

       if (index == 0) then

          ! Linear interpolation or extrapolation
          slope = ( radial_grav(index+1) - radial_grav(index) ) / dr
          mag_grav = radial_grav(index) + slope * xi

       else if (index == n1d-1) then

          ! Linear interpolation or extrapolation
          slope = ( radial_grav(index) - radial_grav(index-1) ) / dr
          mag_grav = radial_grav(index) + slope * xi

       else if (index .gt. n1d-1) then
#ifndef AMREX_USE_CUDA
          if (level .eq. 0) then
             print *,'PUT_RADIAL_GRAV: INDEX TOO BIG ', index, ' > ', n1d-1
             print *,'AT i ', i
             print *,'R / DR IS ', r, dr
             call amrex_error("Error:: Gravity_1d.f90 :: ca_put_radial_grav")
          else
             ! NOTE: we don't do anything to this point if it's outside the
             !       radial grid and level > 0
          end if
#endif
       else

          ! Quadratic interpolation
          ghi = radial_grav(index+1)
          gmd = radial_grav(index  )
          glo = radial_grav(index-1)
          mag_grav = ( ghi -  TWO*gmd + glo)*xi**2/(TWO*dr**2) + &
                     ( ghi       - glo     )*xi   /(TWO*dr   ) + &
                     (-ghi + 26.e0_rt*gmd - glo)/24.e0_rt

          minvar = min(gmd, min(glo,ghi))
          maxvar = max(gmd, max(glo,ghi))
          mag_grav = max(mag_grav,minvar)
          mag_grav = min(mag_grav,maxvar)

       end if

       if (index .le. n1d-1) then

          grav(i,j,k,1) = mag_grav

       end if

    enddo

  end subroutine ca_put_radial_grav



  subroutine ca_put_radial_phi (lo,hi,domlo,domhi,dx,dr,&
       phi,p_lo,p_hi, &
       radial_phi,problo,&
       numpts_1d,fill_interior) bind(C, name="ca_put_radial_phi")

    use amrex_constants_module
    use prob_params_module, only: center

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer , intent(in   ) :: lo(3),hi(3)
    integer , intent(in   ) :: domlo(3),domhi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), value, intent(in   ) :: dr
    real(rt), intent(in   ) :: problo(3)

    integer , value, intent(in   ) :: numpts_1d
    real(rt), intent(in   ) :: radial_phi(0:numpts_1d-1)
    integer , value, intent(in   ) :: fill_interior

    integer , intent(in   ) :: p_lo(3),p_hi(3)
    real(rt), intent(inout) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))

    integer          :: i, j, k, index
    real(rt)         :: r
    real(rt)         :: cen, xi, slope, ph_lo, p_md, ph_hi, minvar, maxvar

    !$gpu

    ! Note that when we interpolate into the ghost cells we use the
    ! location of the edge, not the cell center.

    j = lo(2)
    k = lo(3)

    do i = lo(1), hi(1)
       if (i .gt. domhi(1)) then
          r = problo(1) + (dble(i  )     ) * dx(1) - center(1)
       else if (i .lt. domlo(1)) then
          r = problo(1) + (dble(i+1)     ) * dx(1) - center(1)
       else
          r = problo(1) + (dble(i  )+HALF) * dx(1) - center(1)
       end if

       index = int(r / dr)

#ifndef AMREX_USE_CUDA
       if (index .gt. numpts_1d-1) then
          print *,'PUT_RADIAL_PHI: INDEX TOO BIG ', index, ' > ', numpts_1d-1
          print *,'AT (i) ',i
          print *,'R / DR IS ', r, dr
          call amrex_error("Error:: Gravity_1d.f90 :: ca_put_radial_phi")
       end if
#endif

       if ( (fill_interior .eq. 1) .or. (i .lt. domlo(1) .or. i.gt.domhi(1)) ) then

          cen = (dble(index) + HALF) * dr
          xi = r - cen

          if (index == 0) then

             ! Linear interpolation or extrapolation
             slope = ( radial_phi(index+1) - radial_phi(index) ) / dr
             phi(i,j,k) = radial_phi(index) + slope * xi

          else if (index == numpts_1d-1) then

             ! Linear interpolation or extrapolation
             slope = ( radial_phi(index) - radial_phi(index-1) ) / dr
             phi(i,j,k) = radial_phi(index) + slope * xi

          else

             ! Quadratic interpolation
             ph_hi = radial_phi(index+1)
             p_md = radial_phi(index  )
             ph_lo = radial_phi(index-1)
             phi(i,j,k) = ( ph_hi -  TWO*p_md + ph_lo)*xi**2/(TWO*dr**2) + &
                      ( ph_hi       - ph_lo      )*xi   /(TWO*dr   ) + &
                      (-ph_hi + 26.e0_rt*p_md - ph_lo)/24.e0_rt
             minvar = min(p_md, min(ph_lo,ph_hi))
             maxvar = max(p_md, max(ph_lo,ph_hi))
             phi(i,j,k) = max(phi(i,j,k),minvar)
             phi(i,j,k) = min(phi(i,j,k),maxvar)

          end if

       end if

    enddo

  end subroutine ca_put_radial_phi

end module gravity_1D_module
