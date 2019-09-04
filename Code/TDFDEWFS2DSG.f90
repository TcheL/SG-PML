MODULE InputPara

  IMPLICIT NONE

  PUBLIC
  INTEGER, PARAMETER :: nx = 160, nz = 160                      ! nx: the total number of grid nodes in x-direction; nz: the total number of grid nodes in z-direction.
  INTEGER, PARAMETER :: npmlx = 20, npmlz = 20                  ! npmlx: the total number of grid nodes in top and bottom side of PML absorbing boundary; npmlz: the total number of grid nodes in left and right side of PML absorbing boundary.
  INTEGER, PARAMETER :: sx = 80, sz = 80                        ! sx: the grid node number of source position in x-direction; sz: the grid node number of source position in z-direction.
  INTEGER, PARAMETER :: dx = 5, dz = 5                          ! dx: the grid node interval in x-direction; dz: the grid node interval in z-direction; Unit: m.
  INTEGER, PARAMETER :: nt = 500                                ! the total number of time nodes for wave calculating.
  REAL   , PARAMETER :: dt = 1.0E-3                             ! the time node interval, Unit: s.
  INTEGER, PARAMETER :: nppw = 12                               ! the total node point number per wavelength for dominant frequency of Ricker wavelet source.
  REAL   , PARAMETER :: amp = 1.0E0                             ! the amplitude of source wavelet.
  INTEGER, PARAMETER :: nodr = 3                                ! half of the order number for spatial difference.
  INTEGER, PARAMETER :: irstr = 1                               ! the node ID of starting reciver point.
  INTEGER, PARAMETER :: nrintv = 3                              ! the total node number between each two adjacent recivers.
  INTEGER, PARAMETER :: itstr = 1                               ! the time node ID of the first snapshot.
  INTEGER, PARAMETER :: ntintv = 5                              ! the total time node number between each two followed snapshot.

  REAL    :: src(nt)                                            ! the time series of source wavelet.
  REAL    :: vp(nz, nx), vs(nz, nx), rho(nz, nx)                ! vp: the velocity of P-wave of model, Unit: m/s; vs: the velocity of S-wave of model, Unit: m/s; rho: the density of model, Unit: kg/m^3.   
  INTEGER, PARAMETER :: nrcvr = CEILING(REAL(nx)/nrintv)        ! the total number of all recivers.
  INTEGER :: xrcvr(nrcvr)                                       ! the grid node number in x-direction of reciver position on ground.

  INTEGER, PRIVATE :: i
  
  PRIVATE nppw, amp
  PRIVATE ModelVpRho, SrcWavelet

  CONTAINS
    SUBROUTINE IntlzInputPara()
      xrcvr = [ (irstr + (i - 1)*nrintv, i = 1, nrcvr) ]
      CALL ModelVpRho()
      CALL SrcWavelet()
    END SUBROUTINE IntlzInputPara
    SUBROUTINE ModelVpRho()
      ! here you can reset $vp$ and $rho$ for the model.
      vp = 2000
      vs = 1000
      rho = 1000
    END SUBROUTINE ModelVpRho
    SUBROUTINE SrcWavelet()
      ! here you can reset $src$ for the source wavelet.
      REAL :: f0, t0, pi = 3.1415926
      REAL :: t(nt)
      f0 = MINVAL(vs)/(MIN(dx, dz)*nppw)
      t0 = 1/f0
      t = [ (i*dt, i = 1, nt) ]
      src = amp*(1 - 2*(pi*f0*(t - t0))**2)*EXP( - (pi*f0*(t - t0))**2)
    END SUBROUTINE SrcWavelet

END MODULE InputPara    

MODULE WaveExtrp

  USE InputPara
  IMPLICIT NONE

  PRIVATE
  REAL :: C(nodr)                                               ! the difference coefficients of spatial the $2*nodr$-th order difference approximating.
  INTEGER, PARAMETER :: Nzz = nz + 2*npmlz, Nxx = nx + 2*npmlx  ! Nzz: the total number of grid nodes in z-direction of compute-updating zone including PML layer; Nxx: the total number of grid nodes in x-direction of compute-updating zone including PML layer.
  REAL :: vpp(Nzz, Nxx), vss(Nzz, Nxx), rhoo(Nzz, Nxx)          ! vpp: the velocity of P-wave of the expanded model including PML layer, Unit: m/s; vss: the velocity of S-wave of the expanded model including PML layer, Unit: m/s; rhoo: the density of the expanded model including PML layer, Unit: kg/m^3.
  REAL :: lmdd(Nzz, Nxx), muu(Nzz, Nxx)                         ! lmdd: the lame parameter lambda of elastic wave of the expanded model including PML layer; muu: the lame parameter mu of elastic wave of the expanded model including PML layer.
  REAL :: dpmlz(Nzz, Nxx), dpmlx(Nzz, Nxx)                      ! dpmlz: the PML damping factor in z-direction; dpmlx: the PML damping factor in x-direction.
  REAL :: Coef1(Nzz, Nxx), Coef2(Nzz, Nxx), &
    & Coef3(Nzz, Nxx), Coef4(Nzz, Nxx), &
    & Coef5(Nzz, Nxx), Coef6(Nzz, Nxx), &
    & Coef7(Nzz, Nxx), Coef8(Nzz, Nxx), &
    & Coef9(Nzz, Nxx), Coef0(Nzz, Nxx)                          ! Coef1 ~ Coef0: the coefficients of wavefield time-extrapolating formula.

  INTEGER :: i, j

  REAL, PUBLIC :: P(nz, nx, nt)                                 ! the calculating wavefield component varying with time.

  PUBLIC WaveExec

  CONTAINS
    SUBROUTINE WaveExec()
      CALL CalC()
      CALL ModelExpand()
      CALL CalCoefs()
      CALL CalWave()
    END SUBROUTINE WaveExec
    SUBROUTINE CalC()
      REAL :: rtemp1, rtemp2
      DO i = 1,nodr,1
        rtemp1 = 1.0
        rtemp2 = 1.0
        DO j = 1,nodr,1
          IF(j == i) CYCLE
          rtemp1 = rtemp1*((2*j - 1)**2)
          rtemp2 = rtemp2*ABS((2*i - 1)**2 - (2*j - 1)**2)
        END DO
        C(i) = ( - 1)**(i + 1)*rtemp1/((2*i - 1)*rtemp2)
      END DO
    END SUBROUTINE CalC
    SUBROUTINE ModelExpand()
      vpp = 0.0
      vss = 0.0
      rhoo = 0.0
      vpp(npmlz + 1:npmlz + nz, npmlx + 1:npmlx + nx) = vp
      vss(npmlz + 1:npmlz + nz, npmlx + 1:npmlx + nx) = vs
      rhoo(npmlz + 1:npmlz + nz, npmlx + 1:npmlx + nx) = rho
      DO i = 1,npmlx,1
        vpp(:, i) = vpp(:, npmlx + 1)
        vpp(:, npmlx + nx + i) = vpp(:, npmlx + nx)
        vss(:, i) = vss(:, npmlx + 1)
        vss(:, npmlx + nx + i) = vss(:, npmlx + nx)
        rhoo(:, i) = rhoo(:, npmlx + 1)
        rhoo(:, npmlx + nx + i) = rhoo(:, npmlx + nx)
      END DO
      DO i = 1,npmlz,1
        vpp(i, :) = vpp(npmlz + 1, :)
        vpp(npmlz + nz + i, :) = vpp(npmlz + nz, :)
        vss(i, :) = vss(npmlz + 1, :)
        vss(npmlz + nz + i, :) = vss(npmlz + nz, :)
        rhoo(i, :) = rhoo(npmlz + 1, :)
        rhoo(npmlz + nz + i, :) = rhoo(npmlz + nz, :)
      END DO
      lmdd = rhoo*(vpp**2 - 2*(vss**2))
      muu = rhoo*(vss**2)
    END SUBROUTINE ModelExpand
    SUBROUTINE CalDpml()
      REAL :: dpml0z, dpml0x
      dpml0z = 3*MAXVAL(vs)/dz*(8.0/15 - 3.0/100*npmlz + 1.0/1500*(npmlz**2))
      DO i = 1,npmlz,1
        dpmlz(i, :) = dpml0z*((REAL(npmlz - i + 1)/npmlz)**2)
      END DO
      dpmlz(npmlz + nz + 1:Nzz, :) = dpmlz(npmlz:1:-1, :)
      dpml0x = 3*MAXVAL(vs)/dx*(8.0/15 - 3.0/100*npmlx + 1.0/1500*(npmlx**2))
      DO i = 1,npmlx,1
        dpmlx(:, i) = dpml0x*((REAL(npmlx - i + 1)/npmlx)**2)
      END DO
      dpmlx(:, npmlx + nx + 1:Nxx) = dpmlx(:, npmlx:1:-1)
    END SUBROUTINE CalDpml
    SUBROUTINE CalCoefs()
      CALL CalDpml()
      Coef1 = (2 - dt*dpmlx)/(2 + dt*dpmlx)
      Coef2 = (2 - dt*dpmlz)/(2 + dt*dpmlz)
      Coef3 = (2*dt/(2 + dt*dpmlx))/rhoo/dx
      Coef4 = (2*dt/(2 + dt*dpmlz))/rhoo/dz
      Coef5 = (2*dt/(2 + dt*dpmlx))*(lmdd + 2*muu)/dx
      Coef6 = (2*dt/(2 + dt*dpmlz))*lmdd/dz
      Coef7 = (2*dt/(2 + dt*dpmlx))*lmdd/dx
      Coef8 = (2*dt/(2 + dt*dpmlz))*(lmdd + 2*muu)/dz
      Coef9 = (2*dt/(2 + dt*dpmlx))*muu/dx
      Coef0 = (2*dt/(2 + dt*dpmlz))*muu/dz
    END SUBROUTINE CalCoefs
    SUBROUTINE CalWave()
      INTEGER :: it
      INTEGER, PARAMETER :: Nzzz = Nzz + 2*nodr, Nxxx = Nxx + 2*nodr
      INTEGER :: znds(nz) = [ (nodr + npmlz + i, i = 1,nz,1) ], &
        & xnds(nx) = [ (nodr + npmlx + i, i = 1,nx,1) ]
      INTEGER :: Zznds(Nzz) = [ (nodr + i, i = 1,Nzz,1) ], &
        & Xxnds(Nxx) = [ (nodr + i, i = 1,Nxx,1) ]
      INTEGER :: nsrcz = nodr + npmlz + sz, nsrcx = nodr + npmlx + sx
      REAL    :: vxt(Nzzz, Nxxx) = 0, vxx(Nzzz, Nxxx) = 0, &
        & vxz(Nzzz, Nxxx) = 0, vzt(Nzzz, Nxxx) = 0, &
        & vzx(Nzzz, Nxxx) = 0, vzz(Nzzz, Nxxx) = 0, &
        & txxt(Nzzz, Nxxx) = 0, txxx(Nzzz, Nxxx) = 0, &
        & txxz(Nzzz, Nxxx) = 0, tzzt(Nzzz, Nxxx) = 0, &
        & tzzx(Nzzz, Nxxx) = 0, tzzz(Nzzz, Nxxx) = 0, &
        & txzt(Nzzz, Nxxx) = 0, txzx(Nzzz, Nxxx) = 0, &
        & txzz(Nzzz, Nxxx) = 0, SpcSum(Nzz, Nxx) = 0
      DO it = 1,nt,1
        WRITE(*,"(A,G0)") 'The calculating time node is: it = ',it
        !! load source
        txxx(nsrcz, nsrcx) = txxx(nsrcz, nsrcx) + src(it)/4
        txxz(nsrcz, nsrcx) = txxz(nsrcz, nsrcx) + src(it)/4
        tzzx(nsrcz, nsrcx) = tzzx(nsrcz, nsrcx) + src(it)/4
        tzzz(nsrcz, nsrcx) = tzzz(nsrcz, nsrcx) + src(it)/4
        txxt(:, :) = txxx(:, :) + txxz(:, :)
        tzzt(:, :) = tzzx(:, :) + tzzz(:, :)
        P(:, :, it) = txxt(znds, xnds);
!        P(:, :, it) = tzzt(znds, xnds);
!        P(:, :, it) = txzt(znds, xnds);
        !! calculate $v_x$
        SpcSum = 0
        DO i = 1,nodr,1
          SpcSum = SpcSum + C(i)*(txxt(Zznds, Xxnds + i - 1) - txxt(Zznds, Xxnds - i))
        END DO
        vxx(Zznds, Xxnds) = Coef1*vxx(Zznds, Xxnds) + Coef3*SpcSum
        SpcSum = 0
        DO i = 1,nodr,1
          SpcSum = SpcSum + C(i)*(txzt(Zznds + i - 1, Xxnds) - txzt(Zznds - i, Xxnds))
        END DO
        vxz(Zznds, Xxnds) = Coef2*vxz(Zznds, Xxnds) + Coef4*SpcSum
        vxt(:, :) = vxx(:, :) + vxz(:, :)
!        P(:, :, it) = vxt(znds, xnds)
        !! calculate $v_z$
        SpcSum = 0
        DO i = 1,nodr,1
          SpcSum = SpcSum + C(i)*(txzt(Zznds, Xxnds + i) - txzt(Zznds, Xxnds - i + 1))
        END DO
        vzx(Zznds, Xxnds) = Coef1*vzx(Zznds, Xxnds) + Coef3*SpcSum
        SpcSum = 0
        DO i = 1,nodr,1
          SpcSum = SpcSum + C(i)*(tzzt(Zznds + i, Xxnds) - tzzt(Zznds - i + 1, Xxnds))
        END DO
        vzz(Zznds, Xxnds) = Coef2*vzz(Zznds, Xxnds) + Coef4*SpcSum
        vzt(:, :) = vzx(:, :) + vzz(:, :)
!        P(:, :, it) = vzt(znds, xnds)
        !! calculate $\tau_{xx}$ and $\tau_{zz}$
        SpcSum = 0
        DO i = 1,nodr,1
          SpcSum = SpcSum + C(i)*(vxt(Zznds, Xxnds + i) - vxt(Zznds, Xxnds - i + 1))
        END DO
        txxx(Zznds, Xxnds) = Coef1*txxx(Zznds, Xxnds) + Coef5*SpcSum
        tzzx(Zznds, Xxnds) = Coef1*tzzx(Zznds, Xxnds) + Coef7*SpcSum
        SpcSum = 0
        DO i = 1,nodr,1
          SpcSum = SpcSum + C(i)*(vzt(Zznds + i - 1, Xxnds) - vzt(Zznds - i, Xxnds))
        END DO
        txxz(Zznds, Xxnds) = Coef2*txxz(Zznds, Xxnds) + Coef6*SpcSum
        tzzz(Zznds, Xxnds) = Coef1*tzzz(Zznds, Xxnds) + Coef8*SpcSum
        txxt(:, :) = txxx(:, :) + txxz(:, :)
        tzzt(:, :) = tzzx(:, :) + tzzz(:, :)
        !! calculate $\tau_{xz}$
        SpcSum = 0
        DO i = 1,nodr,1
          SpcSum = SpcSum + C(i)*(vzt(Zznds, Xxnds + i - 1) - vzt(Zznds, Xxnds - i))
        END DO
        txzx(Zznds, Xxnds) = Coef1*txzx(Zznds, Xxnds) + Coef9*SpcSum
        SpcSum = 0
        DO i = 1,nodr,1
          SpcSum = SpcSum + C(i)*(vxt(Zznds + i, Xxnds) - vxt(Zznds - i + 1, Xxnds))
        END DO
        txzz(Zznds, Xxnds) = Coef2*txzz(Zznds, Xxnds) + Coef0*SpcSum
        txzt(:, :) = txzx(:, :) + txzz(:, :)
      END DO
    END SUBROUTINE CalWave 

END MODULE WaveExtrp

!************************  TDFDEWFS2DSG  ***********************
! Time Domain Finite Difference Elastic Wave Field Simulating with 2-Dimension Staggered Grid
! Written by Tche. L. from USTC, 2016.7.
! References:
!   Collino and Tsogka, 2001. Geophysics, Application of the perfectly matched absorbing layer model to the linear elastodynamic problem in anisotropic heterogeneous media.
!   Marcinkovich and Olsen, 2003. Journal of Geophysical Research, On the implementation of perfectly mathced layers in a three-dimensional fourth-order velocity-stress finite difference scheme.
!***************************************************************
PROGRAM TDFDEWFS2DSG

  USE InputPara
  USE WaveExtrp
  IMPLICIT NONE

  CHARACTER(LEN = 128) :: SnapFile = './data/Snapshot_****.dat'       ! the snapshot file name template.
  CHARACTER(LEN = 128) :: SyntFile = './data/SyntRcrd.dat'                       ! the synthetic record file name.
  REAL    :: SyntR(nrcvr, nt)
  INTEGER :: i

  CALL IntlzInputPara()
  CALL WaveExec()
  DO i = itstr,nt,ntintv
    WRITE(SnapFile(21:24),"(I4.4)") i
    CALL Output(TRIM(SnapFile), nz, nx, P(:, :, i))
  END DO
  DO i = 1,nt,1
    SyntR(:, i) = P(1, xrcvr, i)
  END DO
  CALL Output(TRIM(SyntFile), nrcvr, nt, SyntR)

END PROGRAM TDFDEWFS2DSG

SUBROUTINE Output(Outfile, M, N, OutA)
  IMPLICIT NONE
  CHARACTER(LEN = *), INTENT(IN) :: Outfile
  INTEGER, INTENT(IN) :: M, N
  REAL, INTENT(IN)    :: OutA(M, N)
  CHARACTER(LEN = 40) :: FmtStr
  INTEGER :: i, j
  INTEGER :: funit
  WRITE(FmtStr,"('(',G0,'E15.6)')") N
  OPEN(NEWUNIT = funit, FILE = Outfile, STATUS = 'UNKNOWN')
    DO i = 1,M,1
      WRITE(funit, FmtStr) (OutA(i, j), j = 1,N,1)
    END DO
  CLOSE(funit)
END SUBROUTINE Output
