module flatten_module

  implicit none

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine uflaten(lo,hi,p,u,v,w,flatn,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3)

    use mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : iorder, small_pres
    use bl_constants_module

    implicit none

    integer lo(3),hi(3)
    integer q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
    
    double precision p(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision u(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision v(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision w(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision flatn(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)

    integer i, j, k, ishft

    double precision denom, zeta, tst, tmp, ftmp

    ! Local arrays
    double precision, pointer :: dp(:,:,:), z(:,:,:), chi(:,:,:)
    
    ! Knobs for detection of strong shock
    double precision, parameter :: shktst = 0.33d0, zcut1 = 0.75d0, zcut2 = 0.85d0, dzcut = ONE/(zcut2-zcut1)

    call bl_allocate(dp ,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(z  ,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(chi,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)

    !$acc data create(dp, z, chi)

    ! x-direction flattening coef
    !$acc parallel loop collapse(3) present(p, u)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2) 
          !dir$ ivdep
          do i = lo(1)-1,hi(1)+1
             dp(i,j,k) = p(i+1,j,k) - p(i-1,j,k)
             denom = max(small_pres,abs(p(i+2,j,k)-p(i-2,j,k)))
             zeta = abs(dp(i,j,k))/denom
             z(i,j,k) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
             if (u(i-1,j,k)-u(i+1,j,k) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif
             tmp = min(p(i+1,j,k),p(i-1,j,k))
             if ((abs(dp(i,j,k))/tmp).gt.shktst) then
                chi(i,j,k) = tst
             else
                chi(i,j,k) = ZERO
             endif
          enddo
          do i = lo(1),hi(1)
             if(dp(i,j,k).gt.ZERO)then
                ishft = 1
             else
                ishft = -1
             endif
             flatn(i,j,k) = ONE - &
                  max(chi(i-ishft,j,k)*z(i-ishft,j,k),chi(i,j,k)*z(i,j,k))
          enddo
       enddo
    enddo
    !$acc end parallel loop

    ! y-direction flattening coef
    !$acc parallel loop collapse(3) present(p, v)
    do k = lo(3),hi(3)
       do j = lo(2)-1,hi(2)+1
          !dir$ ivdep
          do i = lo(1),hi(1)
             dp(i,j,k) = p(i,j+1,k) - p(i,j-1,k)
             denom = max(small_pres,abs(p(i,j+2,k)-p(i,j-2,k)))
             zeta = abs(dp(i,j,k))/denom
             z(i,j,k) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
             if (v(i,j-1,k)-v(i,j+1,k) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif
             tmp = min(p(i,j+1,k),p(i,j-1,k))
             if ((abs(dp(i,j,k))/tmp).gt.shktst) then
                chi(i,j,k) = tst
             else
                chi(i,j,k) = ZERO
             endif
          enddo
       end do
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if(dp(i,j,k).gt.ZERO)then
                ishft = 1
             else
                ishft = -1
             endif
             ftmp = ONE - &
                  max(chi(i,j-ishft,k)*z(i,j-ishft,k),chi(i,j,k)*z(i,j,k))
             flatn(i,j,k) = min( flatn(i,j,k), ftmp )
          enddo
       enddo
    enddo
    !$acc end parallel loop
 
    ! z-direction flattening coef
    !$acc parallel loop collapse(3) present(p, w)
    do k = lo(3)-1,hi(3)+1
       do j = lo(2),hi(2) 
          !dir$ ivdep
          do i = lo(1),hi(1)
             dp(i,j,k) = p(i,j,k+1) - p(i,j,k-1)
             denom = max(small_pres,abs(p(i,j,k+2)-p(i,j,k-2)))
             zeta = abs(dp(i,j,k))/denom
             z(i,j,k) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
             if (w(i,j,k-1)-w(i,j,k+1) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif
             tmp = min(p(i,j,k+1),p(i,j,k-1))
             if ((abs(dp(i,j,k))/tmp).gt.shktst) then
                chi(i,j,k) = tst
             else
                chi(i,j,k) = ZERO
             endif
          enddo
       enddo
    enddo
    !$acc end parallel loop

    !$acc parallel loop collapse(3)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2) 
          do i = lo(1),hi(1)
             if(dp(i,j,k).gt.ZERO)then
                ishft = 1
             else
                ishft = -1
             endif
             ftmp = ONE - &
                  max(chi(i,j,k-ishft)*z(i,j,k-ishft),chi(i,j,k)*z(i,j,k))
             flatn(i,j,k) = min( flatn(i,j,k), ftmp )
          enddo
       enddo
    enddo
    !$acc end parallel loop

    !$acc end data
    
    call bl_deallocate(dp )
    call bl_deallocate(z  )
    call bl_deallocate(chi)

  end subroutine uflaten

end module flatten_module
