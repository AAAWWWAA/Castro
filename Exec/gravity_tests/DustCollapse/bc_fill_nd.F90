module bc_fill_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  subroutine hypfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="hypfill")

    use amrex_filcc_module, only: amrex_filccn
    use amrex_constants_module, only: HALF
#ifndef AMREX_USE_CUDA
    use castro_error_module, only: castro_error
#endif
    use meth_params_module, only: NVAR,UMX,UMY,UMZ
    use prob_params_module, only: center

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ), value :: time

    integer  :: i, j, k, n
    integer  :: ic,jc,kc
    real(rt) :: x,y,z,r
    real(rt) :: xc,yc,zc,rc
    real(rt) :: mom,momc

    !$gpu

    ! Do this for all the variables, but we will overwrite the momenta below
    call amrex_filccn(lo, hi, adv, adv_lo, adv_hi, NVAR, domlo, domhi, delta, xlo, bc)

#if AMREX_SPACEDIM == 3
#ifndef AMREX_USE_CUDA
    if ( (bc(1,1,1) == EXT_DIR .or. bc(1,2,1) == EXT_DIR) .or.  &
         (bc(2,1,1) == EXT_DIR .or. bc(2,2,1) == EXT_DIR) .or. &
         (bc(3,1,1) == EXT_DIR .or. bc(3,2,1) == EXT_DIR) ) then
       call castro_error("NOT SET UP FOR EXT_DIR BCs IN HYPFILL")
    end if
#endif

    ! XLO
    if ( bc(1,1,1) == FOEXTRAP .and. lo(1) < domlo(1)) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)

             y = (dble(j) + HALF) * delta(2) - center(2)
             z = (dble(k) + HALF) * delta(3) - center(3)

             ic = domlo(1)
             xc = (dble(ic) + HALF) * delta(1) - center(1)
             rc = sqrt(xc**2 + y**2 + z**2)

             momc = sqrt(adv(ic,j,k,UMX)**2 + adv(ic,j,k,UMY)**2 + adv(ic,j,k,UMZ)**2)

             do i = lo(1), hi(1)

                if (i < domlo(1)) then

                   x = (dble(i) + HALF) * delta(1) - center(1)
                   r = sqrt(x**2 + y**2 + z**2)

                   mom  = momc * (rc/r)**2

                   ! Project along the normal
                   adv(i,j,k,UMX) =  sign(1.e0_rt,adv(ic,j,k,UMX)) * mom * (x/r)
                   adv(i,j,k,UMY) =  sign(1.e0_rt,adv(ic,j,k,UMY)) * mom * (y/r)
                   adv(i,j,k,UMZ) =  sign(1.e0_rt,adv(ic,j,k,UMZ)) * mom * (z/r)

                end if

             end do
          end do
       end do
    end if

    ! XHI
    if ( bc(1,2,1) == FOEXTRAP .and. hi(1) > domhi(1)) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)

             y = (dble(j) + HALF) * delta(2) - center(2)
             z = (dble(k) + HALF) * delta(3) - center(3)

             ic = domhi(1)
             xc = (dble(ic) + HALF) * delta(1) - center(1)
             rc = sqrt(xc**2 + y**2 + z**2)

             momc = sqrt(adv(ic,j,k,UMX)**2 + adv(ic,j,k,UMY)**2 + adv(ic,j,k,UMZ)**2)

             do i = lo(1), hi(1)

                if (i > domhi(1)) then

                   x = (dble(i) + HALF) * delta(1) - center(1)
                   r = sqrt(x**2 + y**2 + z**2)

                   mom  = momc * (rc/r)**2

                   ! Project along the normal
                   adv(i,j,k,UMX) =  sign(1.e0_rt,adv(ic,j,k,UMX)) * mom * (x/r)
                   adv(i,j,k,UMY) =  sign(1.e0_rt,adv(ic,j,k,UMY)) * mom * (y/r)
                   adv(i,j,k,UMZ) =  sign(1.e0_rt,adv(ic,j,k,UMZ)) * mom * (z/r)

                end if

             end do
          end do
       end do
    end if

    ! YLO
    if ( bc(2,1,1) == FOEXTRAP .and. lo(2) < domlo(2)) then
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)

             x = (dble(i) + HALF) * delta(1) - center(1)
             z = (dble(k) + HALF) * delta(3) - center(3)

             jc = domlo(2)
             yc = (dble(jc) + HALF) * delta(2) - center(2)
             rc = sqrt(x**2 + yc**2 + z**2)

             momc = sqrt(adv(i,jc,k,UMX)**2 + adv(i,jc,k,UMY)**2 + adv(i,jc,k,UMZ)**2)

             do j = lo(2), hi(2)

                if (j < domlo(2)) then

                   y = (dble(j) + HALF) * delta(2) - center(2)
                   r = sqrt(x**2 + y**2 + z**2)

                   mom  = momc * (rc/r)**2

                   ! Project along the normal
                   adv(i,j,k,UMX) =  sign(1.e0_rt,adv(i,jc,k,UMX)) * mom * (x/r)
                   adv(i,j,k,UMY) =  sign(1.e0_rt,adv(i,jc,k,UMY)) * mom * (y/r)
                   adv(i,j,k,UMZ) =  sign(1.e0_rt,adv(i,jc,k,UMZ)) * mom * (z/r)

                end if

             end do
          end do
       end do
    end if

    ! YHI
    if ( bc(2,2,1) == FOEXTRAP .and. hi(2) > domhi(2)) then
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)

             x = (dble(i) + HALF) * delta(1) - center(1)
             z = (dble(k) + HALF) * delta(3) - center(3)

             jc = domhi(2)
             yc = (dble(jc) + HALF) * delta(2) - center(2)
             rc = sqrt(x**2 + yc**2 + z**2)

             momc = sqrt(adv(i,jc,k,UMX)**2 + adv(i,jc,k,UMY)**2 + adv(i,jc,k,UMZ)**2)

             do j = lo(2), hi(2)

                if (j > domhi(2)) then

                   y = (dble(j) + HALF) * delta(2) - center(2)
                   r = sqrt(x**2 + y**2 + z**2)

                   mom  = momc * (rc/r)**2

                   ! Project along the normal
                   adv(i,j,k,UMX) =  sign(1.e0_rt,adv(i,jc,k,UMX)) * mom * (x/r)
                   adv(i,j,k,UMY) =  sign(1.e0_rt,adv(i,jc,k,UMY)) * mom * (y/r)
                   adv(i,j,k,UMZ) =  sign(1.e0_rt,adv(i,jc,k,UMZ)) * mom * (z/r)

                end if

             end do
          end do
       end do
    end if

    ! ZLO
    if ( bc(3,1,1) == FOEXTRAP .and. lo(3) < domlo(3)) then
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             x = (dble(i) + HALF) * delta(1) - center(1)
             y = (dble(j) + HALF) * delta(2) - center(2)

             kc = domlo(3)
             zc = (dble(kc) + HALF) * delta(3) - center(3)
             rc = sqrt(x**2 + y**2 + zc**2)

             momc = sqrt(adv(i,j,kc,UMX)**2 + adv(i,j,kc,UMY)**2 + adv(i,j,kc,UMZ)**2)

             do k = lo(3), hi(3)

                if (k < domlo(3)) then

                   z = (dble(k) + HALF) * delta(3) - center(3)
                   r = sqrt(x**2 + y**2 + z**2)

                   mom  = momc * (rc/r)**2

                   ! Project along the normal
                   adv(i,j,k,UMX) =  sign(1.e0_rt,adv(i,j,kc,UMX)) * mom * (x/r)
                   adv(i,j,k,UMY) =  sign(1.e0_rt,adv(i,j,kc,UMY)) * mom * (y/r)
                   adv(i,j,k,UMZ) =  sign(1.e0_rt,adv(i,j,kc,UMZ)) * mom * (z/r)

                end if

             end do
          end do
       end do
    end if

    ! ZHI
    if ( bc(3,2,1) == FOEXTRAP .and. hi(3) > domhi(3)) then
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             x = (dble(i) + HALF) * delta(1) - center(1)
             y = (dble(j) + HALF) * delta(2) - center(2)

             kc = domhi(3)
             zc = (dble(kc) + HALF) * delta(3) - center(3)
             rc = sqrt(x**2 + y**2 + zc**2)

             momc = sqrt(adv(i,j,kc,UMX)**2 + adv(i,j,kc,UMY)**2 + adv(i,j,kc,UMZ)**2)

             do k = lo(3), hi(3)

                if (k > domhi(3)) then

                   z = (dble(k) + HALF) * delta(3) - center(3)
                   r = sqrt(x**2 + y**2 + z**2)

                   mom  = momc * (rc/r)**2

                   ! Project along the normal
                   adv(i,j,k,UMX) =  sign(1.e0_rt,adv(i,j,kc,UMX)) * mom * (x/r)
                   adv(i,j,k,UMY) =  sign(1.e0_rt,adv(i,j,kc,UMY)) * mom * (y/r)
                   adv(i,j,k,UMZ) =  sign(1.e0_rt,adv(i,j,kc,UMZ)) * mom * (z/r)

                end if

             end do
          end do
       end do
    end if
#endif

  end subroutine hypfill



  subroutine denfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C,name="denfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
    real(rt), intent(in   ), value :: time

    !$gpu

    call amrex_filccn(lo, hi, adv, adv_lo, adv_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine denfill



#ifdef GRAVITY
  subroutine gravxfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc) bind(C,name="gravxfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))
    real(rt), intent(in   ), value :: time

    !$gpu

    call amrex_filccn(lo, hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravxfill


  subroutine gravyfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc) bind(C,name="gravyfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))
    real(rt), intent(in   ), value :: time

    !$gpu

    call amrex_filccn(lo, hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravyfill


  subroutine gravzfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc) bind(C,name="gravzfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))
    real(rt), intent(in   ), value :: time

    !$gpu

    call amrex_filccn(lo, hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravzfill


  subroutine phigravfill(lo, hi, phi, phi_lo, phi_hi, domlo, domhi, delta, xlo, time, bc) bind(C,name="phigravfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(rt), intent(in   ), value :: time

    !$gpu

    call amrex_filccn(lo, hi, phi, phi_lo, phi_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine phigravfill
#endif

end module bc_fill_module
