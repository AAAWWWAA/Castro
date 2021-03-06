#include "AMReX_LO_BCTYPES.H"
#include "AMReX_ArrayLim.H"

module rad_module

  use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS, UFX, NVAR

  use rad_util_module, only : FLDlambda

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)        , parameter :: tiny = 1.e-50_rt
  real(rt)        , parameter :: BIGKR = 1.e25_rt

contains

subroutine multrs(d, &
                  DIMS(dbox), &
                  DIMS(reg), &
                  r, s) bind(C, name="multrs")
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: DIMDEC(dbox)
  integer :: DIMDEC(reg)
  real(rt)         :: d(DIMV(dbox))
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: s(reg_l2:reg_h2)
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        d(i,j) = d(i,j) * r(i) * s(j)
     enddo
  enddo
end subroutine multrs

subroutine sphc(r, s, &
                DIMS(reg), dx) bind(C, name="sphc")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: s(reg_l2:reg_h2)
  real(rt)         :: dx(2)
  real(rt)         :: h1, h2, d1, d2
  integer :: i, j
  h1 = 0.5e0_rt * dx(1)
  h2 = 0.5e0_rt * dx(2)
  d1 = 1.e0_rt / (3.e0_rt * dx(1))
  d2 = 1.e0_rt / dx(2)
  do i = reg_l1, reg_h1
     r(i) = d1 * ((r(i) + h1)**3 - (r(i) - h1)**3)
  enddo
  do j = reg_l2, reg_h2
     s(j) = d2 * (cos(s(j) - h2) - cos(s(j) + h2))
  enddo
end subroutine sphc

subroutine sphe(r, s, n, &
                DIMS(reg), dx) bind(C, name="sphe")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: s(reg_l2:reg_h2)
  integer :: n
  real(rt)         :: dx(2)
  real(rt)         :: h1, h2, d1, d2
  integer :: i, j
  if (n == 0) then
     do i = reg_l1, reg_h1
        r(i) = r(i)**2
     enddo
     h2 = 0.5e0_rt * dx(2)
     d2 = 1.e0_rt / dx(2)
     do j = reg_l2, reg_h2
        s(j) = d2 * (cos(s(j) - h2) - cos(s(j) + h2))
     enddo
  else
     h1 = 0.5e0_rt * dx(1)
     d1 = 1.e0_rt / (3.e0_rt * dx(1))
     do i = reg_l1, reg_h1
        r(i) = d1 * ((r(i) + h1)**3 - (r(i) - h1)**3)
     enddo
     do j = reg_l2, reg_h2
        s(j) = sin(s(j))
     enddo
  endif
end subroutine sphe

subroutine eddfac(efact, &
                  DIMS(rbox), &
                  DIMS(reg), limiter, n) bind(C, name="eddfac")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: n, limiter
  real(rt)         :: efact(DIMV(rbox))
  integer :: i, j
  real(rt)         :: r, lambda
  if (n == 0) then
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1 + 1
           r = efact(i,j)
           lambda = FLDlambda(r,limiter)
           efact(i,j) = lambda + (lambda * r)**2
        enddo
     enddo
  else
     do j = reg_l2, reg_h2 + 1
        do i = reg_l1, reg_h1
           r = efact(i,j)
           lambda = FLDlambda(r,limiter)
           efact(i,j) = lambda + (lambda * r)**2
        enddo
     enddo
  endif
end subroutine eddfac

subroutine anatw2(test, &
                  DIMS(reg), &
                  temp, p, xf, Tc, dx, xlo, lo) bind(C, name="anatw2")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  real(rt)         :: test(DIMV(reg), 0:1)
  real(rt)         :: temp(DIMV(reg))
  real(rt)         :: p, xf, Tc, dx(2), xlo(2)
  integer :: lo(2)
  integer :: i, j
  real(rt)         :: x, y, r2
  do j = reg_l2, reg_h2
     y = xlo(2) + dx(2) * ((j-lo(2)) + 0.5e0_rt)
     do i = reg_l1, reg_h1
        x  = xlo(1) + dx(1) * ((i-lo(1)) + 0.5e0_rt)
        r2 = x*x + y*y
        test(i,j,0) = Tc * max((1.e0_rt-r2/xf**2), 0.e0_rt)**(1.e0_rt/p)
        test(i,j,1) = temp(i,j) - test(i,j,0)
     enddo
  enddo
end subroutine anatw2

! exch contains temp on input:

subroutine cexch( DIMS(reg), &
     exch, DIMS(xbox), &
     er  , DIMS(ebox), &
     fkp , DIMS(kbox), &
     sigma, c) bind(C, name="cexch")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(xbox)
  integer :: DIMDEC(ebox)
  integer :: DIMDEC(kbox)
  real(rt)         :: exch(DIMV(xbox))
  real(rt)         :: er  (DIMV(ebox))
  real(rt)         :: fkp (DIMV(kbox))
  real(rt)         :: sigma, c
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        exch(i,j) = fkp(i,j) * &
             (4.e0_rt * sigma * exch(i,j)**4 &
             - c * er(i,j))
     enddo
  enddo
end subroutine cexch

subroutine ceupdterm( DIMS(reg), relres, absres, &
     frhoes, DIMS(grd), &
     frhoem, eta, etainv, dfo, dfn, exch, dterm, &
     dt, theta) bind(C, name="ceupdterm")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(grd)
  real(rt)         :: frhoes(DIMV(grd))
  real(rt)         :: frhoem(DIMV(grd))
  real(rt)         :: eta(DIMV(grd))
  real(rt)         :: etainv(DIMV(grd))
  real(rt)         :: dfo(DIMV(grd))
  real(rt)         :: dfn(DIMV(grd))
  real(rt)         :: exch(DIMV(grd))
  real(rt)         :: dterm(DIMV(grd))
  real(rt)         :: dt, theta, relres, absres
  real(rt)         :: tmp, chg, tot
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        chg = 0.e0_rt
        tot = 0.e0_rt
        tmp = eta(i,j) * frhoes(i,j) + &
             etainv(i,j) * &
             (frhoem(i,j) - &
             dt * ((1.e0_rt - theta) * &
             (dfo(i,j) - dfn(i,j)) + &
             exch(i,j))) &
             + dt * dterm(i,j)
        chg = abs(tmp - frhoes(i,j))
        tot = abs(frhoes(i,j))
        frhoes(i,j) = tmp
        absres = max(absres, chg)
        relres = max(relres, chg / (tot + tiny))
     enddo
  enddo
end subroutine ceupdterm

! nonconservative form based on delta B
subroutine nceup(DIMS(reg), relres, absres, &
                 frhoes, DIMS(grd), &
                 frhoem, eta, etainv, &
                 er, DIMS(ebox), &
                 dfo, dfn, temp, fkp, cv, &
                 state, DIMS(sb), &
                 sigma, c, dt, theta) bind(C, name="nceup")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(grd)
  integer :: DIMDEC(sb)
  integer :: DIMDEC(ebox)
  real(rt)         :: frhoes(DIMV(grd))
  real(rt)         :: frhoem(DIMV(grd))
  real(rt)         :: eta(DIMV(grd))
  real(rt)         :: etainv(DIMV(grd))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: dfo(DIMV(grd))
  real(rt)         :: dfn(DIMV(grd))
  real(rt)         :: temp(DIMV(grd))
  real(rt)         :: fkp(DIMV(grd))
  real(rt)         :: cv(DIMV(reg))
  real(rt)         :: state(DIMV(sb), NVAR)
  real(rt)         :: sigma, c, dt, theta, relres, absres
  real(rt)         :: tmp, chg, tot, exch, b, db, dbdt, frhocv
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        chg = 0.e0_rt
        tot = 0.e0_rt
        frhocv = state(i,j, URHO) * cv(i,j)
        dbdt = 16.e0_rt * sigma * temp(i,j)**3
        b = 4.e0_rt * sigma * temp(i,j)**4
        exch = fkp(i,j) * (b - c * er(i,j))
        tmp = eta(i,j) * frhoes(i,j) + etainv(i,j) * &
             (frhoem(i,j) - &
             dt * ((1.e0_rt - theta) * &
             (dfo(i,j) - dfn(i,j)) + &
             exch))
#if 1
        if (frhocv > tiny .AND. tmp > frhoes(i,j)) then
           db = (tmp - frhoes(i,j)) * dbdt / frhocv
           if (b + db <= 0.e0_rt) then
              print *, i, j, b, db, b+db
           endif
           tmp = ((b + db) / (4.e0_rt * sigma))**0.25e0_rt
           tmp = frhoes(i,j) + frhocv * (tmp - temp(i,j))
        endif
#endif
        chg = abs(tmp - frhoes(i,j))
        tot = abs(frhoes(i,j))
        frhoes(i,j) = tmp
        absres = max(absres, chg)
        relres = max(relres, chg / (tot + tiny))
     enddo
  enddo
end subroutine nceup

subroutine cetot(DIMS(reg), &
                 state, DIMS(sb), &
                 frhoe, DIMS(fb)) bind(C, name="cetot")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(sb)
  integer :: DIMDEC(fb)
  real(rt)         :: state(DIMV(sb), NVAR)
  real(rt)         :: frhoe(DIMV(fb))
  real(rt)         :: kin
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        !            kin = 0.5e0_rt * (state(i,j,XMOM)   ** 2 +
        !     @                     state(i,j,XMOM+1) ** 2) /
        !     @                    state(i,j,DEN)
        kin = state(i,j, UEDEN) - state(i,j, UEINT)
        state(i,j, UEINT) = frhoe(i,j)
        state(i,j, UEDEN) = frhoe(i,j) + kin
     enddo
  enddo
end subroutine cetot

! *********************************
! ** BEGIN MGFLD routines        **
! *********************************

subroutine lacoefmgfld(a, &
                       DIMS(abox), &
                       DIMS(reg), &
                       kappa, &
                       DIMS(kbox), &
                       r, s, &
                       dt, c) bind(C, name="lacoefmgfld")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(abox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)

  real(rt)         :: a(DIMV(abox))
  real(rt)         :: kappa(DIMV(kbox))
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: s(reg_l2:reg_h2)
  real(rt)         :: dt, c

  integer :: i, j

  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1

        a(i,j) = c*kappa(i,j) + 1.e0_rt/dt
        a(i,j) = r(i) * s(j) * a(i,j)

     enddo
  enddo
end subroutine lacoefmgfld

! *********************************
! ** END MGFLD routines          **
! *********************************

subroutine rfface(fine, &
                  DIMS(fbox), &
                  crse, &
                  DIMS(cbox), &
                  idim, irat) bind(C, name="rfface")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(fbox)
  integer :: DIMDEC(cbox)
  real(rt)         :: fine(DIMV(fbox))
  real(rt)         :: crse(DIMV(cbox))
  integer :: idim, irat(0:1)
  integer :: i, j
  if (idim == 0) then
     do j = fbox_l2, fbox_h2
        fine(fbox_l1,j) = crse(cbox_l1, j/irat(1)) / irat(1)
     enddo
  else
     do i = fbox_l1, fbox_h1
        fine(i,fbox_l2) = crse(i/irat(0), cbox_l2) / irat(0)
     enddo
  endif
end subroutine rfface


subroutine bextrp(f, fbox_l1, fbox_l2, fbox_h1, fbox_h2, &
                  reg_l1, reg_l2, reg_h1, reg_h2) bind(C, name="bextrp")

  use amrex_fort_module, only : rt => amrex_real
  integer :: fbox_l1, fbox_l2, fbox_h1, fbox_h2
  integer ::  reg_l1,  reg_l2,  reg_h1,  reg_h2
  real(rt)         :: f(fbox_l1:fbox_h1,fbox_l2:fbox_h2)
  integer :: i, j

  !     i direction first:
  do j = reg_l2, reg_h2
     i = reg_l1
     f(i-1,j) = 2.e0_rt * f(i,j) - f(i+1,j)
     i = reg_h1
     f(i+1,j) = 2.e0_rt * f(i,j) - f(i-1,j)
  enddo

  !     j direction second, including corners:
  do i = reg_l1 - 1, reg_h1 + 1
     j = reg_l2
     f(i,j-1) = 2.e0_rt * f(i,j) - f(i,j+1)
     j = reg_h2
     f(i,j+1) = 2.e0_rt * f(i,j) - f(i,j-1)
  enddo

  !  corner results are the same whichever direction we extrapolate first
end subroutine bextrp


subroutine lbcoefna(bcoef, &
                    bcgrp, bboxl0, bboxl1, bboxh0, bboxh1, &
                    reg_l1, reg_l2, reg_h1, reg_h2, &
                    spec, sboxl0, sboxl1, sboxh0, sboxh1, &
                    idim) bind(C, name="lbcoefna")

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: idim
  integer ::  reg_l1,  reg_l2,  reg_h1,  reg_h2
  integer :: bboxl0, bboxl1, bboxh0, bboxh1
  integer :: sboxl0, sboxl1, sboxh0, sboxh1
  real(rt)         :: bcoef(bboxl0:bboxh0,bboxl1:bboxh1)
  real(rt)         :: bcgrp(bboxl0:bboxh0,bboxl1:bboxh1)
  real(rt)         :: spec(sboxl0:sboxh0,sboxl1:sboxh1)
  integer :: i, j
  if (idim == 0) then
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           bcoef(i,j) = bcoef(i,j) &
                + 0.5e0_rt * (spec(i-1,j) + spec(i,j)) * bcgrp(i,j)
        enddo
     enddo
  else
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           bcoef(i,j) = bcoef(i,j) &
                + 0.5e0_rt * (spec(i,j-1) + spec(i,j)) * bcgrp(i,j)
        enddo
     enddo
  endif

end subroutine lbcoefna


subroutine ljupna(jnew, jboxl0, jboxl1, jboxh0, jboxh1, &
                  reg_l1, reg_l2, reg_h1, reg_h2, &
                  spec, sboxl0, sboxl1, sboxh0, sboxh1, &
                  accel, aboxl0, aboxl1, aboxh0, aboxh1, &
                  nTotal) bind(C, name="ljupna")

  use amrex_fort_module, only : rt => amrex_real
  integer :: nTotal
  integer ::  reg_l1,  reg_l2,  reg_h1, reg_h2
  integer :: jboxl0, jboxl1, jboxh0, jboxh1
  integer :: sboxl0, sboxl1, sboxh0, sboxh1
  integer :: aboxl0, aboxl1, aboxh0, aboxh1
  real(rt)         :: jnew(jboxl0:jboxh0,jboxl1:jboxh1,0:nTotal-1)
  real(rt)         :: spec(sboxl0:sboxh0,sboxl1:sboxh1,0:nTotal-1)
  real(rt)         :: accel(aboxl0:aboxh0,aboxl1:aboxh1)

  integer :: i, j, n

  do n = 0, nTotal - 1
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           jnew(i,j,n) = jnew(i,j,n) + spec(i,j,n) * accel(i,j)
        enddo
     enddo
  enddo

end subroutine ljupna

end module rad_module
