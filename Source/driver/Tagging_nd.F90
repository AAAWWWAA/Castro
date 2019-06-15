module tagging_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), allocatable ::    denerr,   dengrad, dengrad_rel
  real(rt), allocatable ::    enterr,   entgrad, entgrad_rel
  real(rt), allocatable ::    velerr,   velgrad, velgrad_rel
  real(rt), allocatable ::   temperr,  tempgrad, tempgrad_rel
  real(rt), allocatable ::  presserr, pressgrad, pressgrad_rel
  real(rt), allocatable ::    raderr,   radgrad, radgrad_rel
  real(rt), allocatable ::   enucerr

  integer, allocatable ::  max_denerr_lev,   max_dengrad_lev, max_dengrad_rel_lev
  integer, allocatable ::  max_velerr_lev,   max_velgrad_lev, max_velgrad_rel_lev
  integer, allocatable ::  max_temperr_lev,  max_tempgrad_lev, max_tempgrad_rel_lev
  integer, allocatable ::  max_presserr_lev, max_pressgrad_lev, max_pressgrad_rel_lev
  integer, allocatable ::  max_raderr_lev,   max_radgrad_lev, max_radgrad_rel_lev
  integer, allocatable ::  max_enucerr_lev

  ! limit the zone size based on how much the burning can change the
  ! internal energy of a zone. The zone size on the finest level must
  ! be smaller than dxnuc * c_s * (e/ \dot{e}) where c_s is the sound
  ! speed.  This ensures that the sound-crossing time is smaller than
  ! the nuclear energy injection timescale.
  real(rt), allocatable :: dxnuc_min

  ! Disable limiting based on dxnuc above this threshold. This allows
  !  zones that have already ignited or are about to ignite to be
  !  de-refined.
  real(rt), allocatable :: dxnuc_max

  ! Disable limiting based on dxnuc above this AMR level.
  integer, allocatable :: max_dxnuc_lev

  public

#ifdef AMREX_USE_CUDA
  attributes(managed) ::    denerr,   dengrad, dengrad_rel
  attributes(managed) ::    velerr,   velgrad, velgrad_rel
  attributes(managed) ::   temperr,  tempgrad, tempgrad_rel
  attributes(managed) ::  presserr, pressgrad, pressgrad_rel
  attributes(managed) ::    raderr,   radgrad, radgrad_rel
  attributes(managed) ::   enucerr

  attributes(managed) ::  max_denerr_lev,   max_dengrad_lev, max_dengrad_rel_lev
  attributes(managed) ::  max_velerr_lev,   max_velgrad_lev, max_velgrad_rel_lev
  attributes(managed) ::  max_temperr_lev,  max_tempgrad_lev, max_tempgrad_rel_lev
  attributes(managed) ::  max_presserr_lev, max_pressgrad_lev, max_pressgrad_rel_lev
  attributes(managed) ::  max_raderr_lev,   max_radgrad_lev, max_radgrad_rel_lev
  attributes(managed) ::  max_enucerr_lev

  attributes(managed) :: dxnuc_min
  attributes(managed) :: dxnuc_max
  attributes(managed) :: max_dxnuc_lev
#endif

contains

  subroutine ca_denerror(lo, hi, &
                         tag, taglo, taghi, &
                         state, slo, shi, &
                         dx, problo, &
                         set, clear, time, level) &
                         bind(C, name="ca_denerror")
    !
    ! This routine will tag high error cells based on the density
    !

    use amrex_constants_module, only: ONE
    use meth_params_module, only: NVAR, URHO, UTEMP, UMX
    use prob_params_module, only: dg

    implicit none

    integer,    intent(in   ) :: lo(3), hi(3)
    integer,    intent(in   ) :: taglo(3), taghi(3)
    integer,    intent(in   ) :: slo(3), shi(3)
    integer(1), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt),   intent(in   ) :: state(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),NVAR)
    real(rt),   intent(in   ) :: dx(3), problo(3)
    integer(1), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: level
    real(rt),   intent(in   ), value :: time

    real(rt) :: ax, ay, az
    real(rt) :: rhoInv
    integer  :: i, j, k, n

    !$gpu

    ! Tag on regions of high density
    if (level .lt. max_denerr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (state(i,j,k,URHO) .ge. denerr) then
                   tag(i,j,k) = set
                end if
             end do
          end do
       end do
    end if

    ! Tag on regions of high density gradient
    if (level .lt. max_dengrad_lev .or. level .lt. max_dengrad_rel_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(state(i+1*dg(1),j,k,URHO) - state(i,j,k,URHO))
                ay = ABS(state(i,j+1*dg(2),k,URHO) - state(i,j,k,URHO))
                az = ABS(state(i,j,k+1*dg(3),URHO) - state(i,j,k,URHO))
                ax = MAX(ax,ABS(state(i,j,k,URHO) - state(i-1*dg(1),j,k,URHO)))
                ay = MAX(ay,ABS(state(i,j,k,URHO) - state(i,j-1*dg(2),k,URHO)))
                az = MAX(az,ABS(state(i,j,k,URHO) - state(i,j,k-1*dg(3),URHO)))
                if (MAX(ax,ay,az) .ge. dengrad .or. MAX(ax,ay,az) .ge. ABS(dengrad_rel * state(i,j,k,URHO))) then
                   tag(i,j,k) = set
                end if
             end do
          end do
       end do
    end if

    ! Tag on regions of high temperature
    if (level .lt. max_temperr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (state(i,j,k,UTEMP) .ge. temperr) then
                   tag(i,j,k) = set
                end if
             end do
          end do
       end do
    end if

    ! Tag on regions of high temperature gradient
    if (level .lt. max_tempgrad_lev .or. level .lt. max_tempgrad_rel_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(state(i+1*dg(1),j,k,UTEMP) - state(i,j,k,UTEMP))
                ay = ABS(state(i,j+1*dg(2),k,UTEMP) - state(i,j,k,UTEMP))
                az = ABS(state(i,j,k+1*dg(3),UTEMP) - state(i,j,k,UTEMP))
                ax = MAX(ax,ABS(state(i,j,k,UTEMP) - state(i-1*dg(1),j,k,UTEMP)))
                ay = MAX(ay,ABS(state(i,j,k,UTEMP) - state(i,j-1*dg(2),k,UTEMP)))
                az = MAX(az,ABS(state(i,j,k,UTEMP) - state(i,j,k-1*dg(3),UTEMP)))
                if (MAX(ax,ay,az) .ge. tempgrad .or. MAX(ax,ay,az) .ge. ABS(tempgrad_rel * state(i,j,k,UTEMP))) then
                   tag(i,j,k) = set
                end if
             end do
          end do
       end do
    end if

    ! Tag on regions of high velocity
    if (level .lt. max_velerr_lev) then
       do n = UMX, UMX + AMREX_SPACEDIM - 1
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if (state(i,j,k,n) .ge. velerr * state(i,j,k,URHO)) then
                      tag(i,j,k) = set
                   end if
                end do
             end do
          end do
       end do
    end if

    ! Tag on regions of high velocity gradient
    if (level .lt. max_velgrad_lev .or. level .lt. max_velgrad_rel_lev) then
       do n = UMX, UMX + AMREX_SPACEDIM - 1
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   rhoInv = ONE / state(i,j,k,URHO)

                   ax = ABS(state(i+1*dg(1),j,k,n) / state(i+1*dg(1),j,k,URHO) - state(i,j,k,n) * rhoInv)
                   ay = ABS(state(i,j+1*dg(2),k,n) / state(i,j+1*dg(2),k,URHO) - state(i,j,k,n) * rhoInv)
                   az = ABS(state(i,j,k+1*dg(3),n) / state(i,j,k+1*dg(3),URHO) - state(i,j,k,n) * rhoInv)
                   ax = MAX(ax,ABS(state(i,j,k,n) * rhoInv - state(i-1*dg(1),j,k,n) / state(i-1*dg(1),j,k,URHO)))
                   ay = MAX(ay,ABS(state(i,j,k,n) * rhoInv - state(i,j-1*dg(2),k,n) / state(i,j-1*dg(2),k,URHO)))
                   az = MAX(az,ABS(state(i,j,k,n) * rhoInv - state(i,j,k-1*dg(3),n) / state(i,j,k-1*dg(3),URHO)))
                   if (MAX(ax,ay,az) .ge. velgrad .or. MAX(ax,ay,az) .ge. ABS(velgrad_rel * state(i,j,k,n) * rhoInv)) then
                      tag(i,j,k) = set
                   end if
                end do
             end do
          end do
       end do
    end if

  end subroutine ca_denerror



  subroutine ca_presserror(lo, hi, &
                           tag, taglo, taghi, &
                           press, presslo, presshi, np, &
                           delta, problo, &
                           set, clear, time, level) &
                           bind(C, name="ca_presserror")
   !
   ! This routine will tag high error cells based on the pressure
   !

    use prob_params_module, only: dg
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,    intent(in   ) :: lo(3), hi(3)
    integer,    intent(in   ) :: taglo(3), taghi(3)
    integer,    intent(in   ) :: presslo(3), presshi(3)
    integer(1), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt),   intent(in   ) :: press(presslo(1):presshi(1),presslo(2):presshi(2),presslo(3):presshi(3),np)
    real(rt),   intent(in   ) :: delta(3), problo(3)
    integer(1), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: np, level
    real(rt),   intent(in   ), value :: time

    real(rt) :: ax, ay, az
    integer  :: i, j, k

    !$gpu

    !     Tag on regions of high pressure
    if (level .lt. max_presserr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (press(i,j,k,1) .ge. presserr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high pressure gradient
    if (level .lt. max_pressgrad_lev .or. level .lt. max_pressgrad_rel_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(press(i+1*dg(1),j,k,1) - press(i,j,k,1))
                ay = ABS(press(i,j+1*dg(2),k,1) - press(i,j,k,1))
                az = ABS(press(i,j,k+1*dg(3),1) - press(i,j,k,1))
                ax = MAX(ax,ABS(press(i,j,k,1) - press(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(press(i,j,k,1) - press(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(press(i,j,k,1) - press(i,j,k-1*dg(3),1)))
                if (MAX(ax,ay,az) .ge. pressgrad .or. MAX(ax,ay,az) .ge. ABS(pressgrad_rel * press(i,j,k,1))) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine ca_presserror



  subroutine ca_raderror(lo, hi, &
                         tag, taglo, taghi, &
                         rad, radlo, radhi, nr, &
                         delta, problo, &
                         set, clear, time, level) &
                         bind(C, name="ca_raderror")
    !
    ! This routine will tag high error cells based on the radiation
    !

    use prob_params_module, only: dg
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,    intent(in   ) :: lo(3), hi(3)
    integer,    intent(in   ) :: taglo(3), taghi(3)
    integer,    intent(in   ) :: radlo(3), radhi(3)
    integer(1), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt),   intent(in   ) :: rad(radlo(1):radhi(1),radlo(2):radhi(2),radlo(3):radhi(3),nr)
    real(rt),   intent(in   ) :: delta(3), problo(3)
    integer(1), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: nr, level
    real(rt),   intent(in   ), value :: time

    real(rt) :: ax, ay, az
    integer  :: i, j, k

    !$gpu

    !     Tag on regions of high radiation
    if (level .lt. max_raderr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (rad(i,j,k,1) .ge. raderr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high radiation gradient
    if (level .lt. max_radgrad_lev .or. level .lt. max_radgrad_rel_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(rad(i+1*dg(1),j,k,1) - rad(i,j,k,1))
                ay = ABS(rad(i,j+1*dg(2),k,1) - rad(i,j,k,1))
                az = ABS(rad(i,j,k+1*dg(3),1) - rad(i,j,k,1))
                ax = MAX(ax,ABS(rad(i,j,k,1) - rad(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(rad(i,j,k,1) - rad(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(rad(i,j,k,1) - rad(i,j,k-1*dg(3),1)))
                if (MAX(ax,ay,az) .ge. radgrad .or. MAX(ax,ay,az) .ge. ABS(radgrad_rel * rad(i,j,k,1))) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine ca_raderror



#ifdef REACTIONS
  subroutine ca_nucerror(lo, hi, &
                         tag, taglo, taghi, &
                         t, tlo, thi, nr, &
                         delta, problo, &
                         set, clear, time, level) &
                         bind(C, name="ca_nucerror")
    !
    ! This routine will tag cells based on the sound crossing time
    ! relative to the nuclear energy injection timescale.
    !

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,    intent(in   ) :: lo(3), hi(3)
    integer,    intent(in   ) :: taglo(3), taghi(3)
    integer,    intent(in   ) :: tlo(3), thi(3)
    integer(1), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt),   intent(in   ) :: t(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),nr) ! t_sound / t_e
    real(rt),   intent(in   ) :: delta(3), problo(3)
    integer(1), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: nr, level
    real(rt),   intent(in   ), value :: time

    integer :: i, j, k

    !$gpu

    ! Disable if we're not utilizing this tagging

    if (dxnuc_min > 1.e199_rt) return

    if (level .lt. max_dxnuc_lev) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (t(i,j,k,1) > dxnuc_min .and. t(i,j,k,1) < dxnuc_max) then

                   tag(i,j,k) = set

                endif

             enddo
          enddo
       enddo

    end if

  end subroutine ca_nucerror



  subroutine ca_enucerror(lo, hi, &
                          tag, taglo, taghi, &
                          reactions, r_lo, r_hi, &
                          delta, problo, &
                          set, clear, time, level) &
                          bind(C, name="ca_enucerror")
    !
    ! This routine will tag high error cells based on the nuclear
    ! energy generation rate
    !

    use prob_params_module, only: dg
    use amrex_fort_module, only: rt => amrex_real
    use network, only: nspec

    implicit none

    integer,    intent(in   ) :: lo(3), hi(3)
    integer,    intent(in   ) :: taglo(3), taghi(3)
    integer,    intent(in   ) :: r_lo(3), r_hi(3)
    integer(1), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt),   intent(in   ) :: reactions(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nspec+2)
    real(rt),   intent(in   ) :: delta(3), problo(3)
    integer(1), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: nd, level
    real(rt),   intent(in   ), value :: time

    integer :: i, j, k

    !$gpu

    ! Tag on regions of high nuclear energy generation rate
    if (level .lt. max_enucerr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (reactions(i,j,k,nspec+1) .ge. enucerr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine ca_enucerror
#endif



  ! Routines for retrieving the maximum tagging level.

  subroutine get_max_denerr_lev(lev) bind(c, name='get_max_denerr_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_denerr_lev

  end subroutine get_max_denerr_lev



  subroutine get_max_dengrad_lev(lev) bind(c, name='get_max_dengrad_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_dengrad_lev

  end subroutine get_max_dengrad_lev



  subroutine get_max_dengrad_rel_lev(lev) bind(c, name='get_max_dengrad_rel_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_dengrad_rel_lev

  end subroutine get_max_dengrad_rel_lev



  subroutine get_max_velerr_lev(lev) bind(c, name='get_max_velerr_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_velerr_lev

  end subroutine get_max_velerr_lev



  subroutine get_max_velgrad_lev(lev) bind(c, name='get_max_velgrad_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_velgrad_lev

  end subroutine get_max_velgrad_lev



  subroutine get_max_velgrad_rel_lev(lev) bind(c, name='get_max_velgrad_rel_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_velgrad_rel_lev

  end subroutine get_max_velgrad_rel_lev



  subroutine get_max_temperr_lev(lev) bind(c, name='get_max_temperr_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_temperr_lev

  end subroutine get_max_temperr_lev



  subroutine get_max_tempgrad_lev(lev) bind(c, name='get_max_tempgrad_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_tempgrad_lev

  end subroutine get_max_tempgrad_lev



  subroutine get_max_tempgrad_rel_lev(lev) bind(c, name='get_max_tempgrad_rel_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_tempgrad_rel_lev

  end subroutine get_max_tempgrad_rel_lev



  subroutine get_max_presserr_lev(lev) bind(c, name='get_max_presserr_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_presserr_lev

  end subroutine get_max_presserr_lev



  subroutine get_max_pressgrad_lev(lev) bind(c, name='get_max_pressgrad_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_pressgrad_lev

  end subroutine get_max_pressgrad_lev



  subroutine get_max_pressgrad_rel_lev(lev) bind(c, name='get_max_pressgrad_rel_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_pressgrad_rel_lev

  end subroutine get_max_pressgrad_rel_lev



  subroutine get_max_raderr_lev(lev) bind(c, name='get_max_raderr_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_raderr_lev

  end subroutine get_max_raderr_lev



  subroutine get_max_radgrad_lev(lev) bind(c, name='get_max_radgrad_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_radgrad_lev

  end subroutine get_max_radgrad_lev



  subroutine get_max_radgrad_rel_lev(lev) bind(c, name='get_max_radgrad_rel_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_radgrad_rel_lev

  end subroutine get_max_radgrad_rel_lev

end module tagging_module
