MODULE InputPara

  IMPLICIT NONE

  PUBLIC
  INTEGER, PARAMETER :: nx = 101, nz = 101                      ! nx: the total number of grid nodes in x-direction; nz: the total number of grid nodes in z-direction.
  INTEGER, PARAMETER :: npmlx = 20, npmlz = 20                  ! npmlx: the total number of grid nodes in top and bottom side of PML absorbing boundary; npmlz: the total number of grid nodes in left and right side of PML absorbing boundary.
  INTEGER, PARAMETER :: sx = 50, sz = 50                        ! sx: the grid node number of source position in x-direction; sz: the grid node number of source position in z-direction.
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
  REAL    :: vp(nz, nx), rho(nz, nx)                            ! vp: the velocity of acoustic wave of model, Unit: m/s; rho: the density of model, Unit: kg/m^3.   
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
      vp(nz/3:nz, nx/2:nx) = 1000
      rho = 1000
      rho(nz/3:nz, nx/2:nx) = 500
    END SUBROUTINE ModelVpRho
    SUBROUTINE SrcWavelet()
      ! here you can reset $src$ for the source wavelet.
      REAL :: f0, t0, pi = 3.1415926
      REAL :: t(nt)
      f0 = MINVAL(vp)/(MIN(dx, dz)*nppw)
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
  REAL :: vpp(Nzz, Nxx), rhoo(Nzz, Nxx)                         ! vpp: the velocity of the expanded model including PML layer, Unit: m/s; rhoo: the density of the expanded model including PML layer, Unit: kg/m^3.
  REAL :: dpmlz(Nzz, Nxx), dpmlx(Nzz, Nxx)                      ! dpmlz: the PML damping factor in z-direction; dpmlx: the PML damping factor in x-direction.
  REAL :: Coef1(Nzz, Nxx), Coef2(Nzz, Nxx), &
    & Coef3(Nzz, Nxx), Coef4(Nzz, Nxx), &
    & Coef5(Nzz, Nxx), Coef6(Nzz, Nxx)                          ! Coef1 ~ Coef6: the coefficients of wavefield time-extrapolating formula.

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
      rhoo = 0.0
      vpp(npmlz + 1:npmlz + nz, npmlx + 1:npmlx + nx) = vp
      rhoo(npmlz + 1:npmlz + nz, npmlx + 1:npmlx + nx) = rho
      DO i = 1,npmlx,1
        vpp(:, i) = vpp(:, npmlx + 1)
        vpp(:, npmlx + nx + i) = vpp(:, npmlx + nx)
        rhoo(:, i) = rhoo(:, npmlx + 1)
        rhoo(:, npmlx + nx + i) = rhoo(:, npmlx + nx)
      END DO
      DO i = 1,npmlz,1
        vpp(i, :) = vpp(npmlz + 1, :)
        vpp(npmlz + nz + i, :) = vpp(npmlz + nz, :)
        rhoo(i, :) = rhoo(npmlz + 1, :)
        rhoo(npmlz + nz + i, :) = rhoo(npmlz + nz, :)
      END DO
    END SUBROUTINE ModelExpand
    SUBROUTINE CalDpml()
      REAL :: dpml0z, dpml0x
      dpml0z = 3*MAXVAL(vp)/dz*(8.0/15 - 3.0/100*npmlz + 1.0/1500*(npmlz**2))
      DO i = 1,npmlz,1
        dpmlz(i, :) = dpml0z*((REAL(npmlz - i + 1)/npmlz)**2)
      END DO
      dpmlz(npmlz + nz + 1:Nzz, :) = dpmlz(npmlz:1:-1, :)
      dpml0x = 3*MAXVAL(vp)/dx*(8.0/15 - 3.0/100*npmlx + 1.0/1500*(npmlx**2))
      DO i = 1,npmlx,1
        dpmlx(:, i) = dpml0x*((REAL(npmlx - i + 1)/npmlx)**2)
      END DO
      dpmlx(:, npmlx + nx + 1:Nxx) = dpmlx(:, npmlx:1:-1)
    END SUBROUTINE CalDpml
    SUBROUTINE CalCoefs()
      CALL CalDpml()
      Coef1 = (2 - dt*dpmlx)/(2 + dt*dpmlx)
      Coef2 = (2 - dt*dpmlz)/(2 + dt*dpmlz)
      Coef3 = (2*dt/(2 + dt*dpmlx))/(rhoo*dx)
      Coef4 = (2*dt/(2 + dt*dpmlz))/(rhoo*dz)
      Coef5 = (2*dt/(2 + dt*dpmlx))*(rhoo*(vpp**2)/dx)
      Coef6 = (2*dt/(2 + dt*dpmlz))*(rhoo*(vpp**2)/dz)
    END SUBROUTINE CalCoefs
    SUBROUTINE CalWave()
      INTEGER :: it
      INTEGER, PARAMETER :: Nzzz = Nzz + 2*nodr, Nxxx = Nxx + 2*nodr
      INTEGER :: znds(nz) = [ (nodr + npmlz + i, i = 1,nz,1) ], &
        & xnds(nx) = [ (nodr + npmlx + i, i = 1,nx,1) ]
      INTEGER :: Zznds(Nzz) = [ (nodr + i, i = 1,Nzz,1) ], &
        & Xxnds(Nxx) = [ (nodr + i, i = 1,Nxx,1) ]
      INTEGER :: nsrcz = nodr + npmlz + sz, nsrcx = nodr + npmlx + sx
      REAL    :: Pt(Nzzz, Nxxx) = 0, &
        & Pz(Nzzz, Nxxx) = 0, Px(Nzzz, Nxxx) = 0, &
        & vz(Nzzz, Nxxx) = 0, vx(Nzzz, Nxxx) = 0
      REAL    :: SpcSum(Nzz, Nxx)
      DO it = 1,nt,1
        WRITE(*,"(A,G0)") 'The calculating time node is: it = ',it
        Px(nsrcz, nsrcx) = Px(nsrcz, nsrcx) + src(it)/2
        Pz(nsrcz, nsrcx) = Pz(nsrcz, nsrcx) + src(it)/2
        Pt(:, :) = Px(:, :) + Pz(:, :)
        P(:, :, it) = Pt(znds, xnds)
        SpcSum = 0
        DO i = 1,nodr,1
          SpcSum = SpcSum + C(i)*(Pt(Zznds, Xxnds + i) - Pt(Zznds, Xxnds + 1 - i))
        END DO
        vx(Zznds, Xxnds) = Coef1*vx(Zznds, Xxnds) - Coef3*SpcSum
        SpcSum = 0
        DO i = 1,nodr,1
          SpcSum = SpcSum + C(i)*(Pt(Zznds + i, Xxnds) - Pt(Zznds + 1 - i, Xxnds))
        END DO
        vz(Zznds, Xxnds) = Coef2*vz(Zznds, Xxnds) - Coef4*SpcSum
        SpcSum = 0
        DO i = 1,nodr,1
          SpcSum = SpcSum + C(i)*(vx(Zznds, Xxnds - 1 + i) - vx(Zznds, Xxnds - i))
        END DO
        Px(Zznds, Xxnds) = Coef1*Px(Zznds, Xxnds) - Coef5*SpcSum
        SpcSum = 0
        DO i = 1,nodr,1
          SpcSum = SpcSum + C(i)*(vz(Zznds - 1 + i, Xxnds) - vz(Zznds - i, Xxnds))
        END DO
        Pz(Zznds, Xxnds) = Coef2*Pz(Zznds, Xxnds) - Coef6*SpcSum
      END DO
    END SUBROUTINE CalWave

END MODULE WaveExtrp

!************************  TDFDAWFS2DSG  ***********************
! Time Domain Finite Difference Acoustic Wave Field Simulating with 2-Dimension Staggered Grid
! Written by Tche. L. from USTC, 2016.7.
! References:
!   Collino and Tsogka, 2001. Geophysics, Application of the perfectly matched absorbing layer model to the linear elastodynamic problem in anisotropic heterogeneous media.
!   Marcinkovich and Olsen, 2003. Journal of Geophysical Research, On the implementation of perfectly mathced layers in a three-dimensional fourth-order velocity-stress finite difference scheme.
!***************************************************************
PROGRAM TDFDAWFS2DSG

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

END PROGRAM TDFDAWFS2DSG

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
