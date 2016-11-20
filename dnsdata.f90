!============================================!
!                                            !
!    Data Structures, Definitions and I/O    !
!                  for the                   !
!     Direct Numerical Simulation (DNS)      !
!        of a turbulent channel flow         !
!                                            !
!============================================!
!
! Author: Dr.-Ing. Davide Gatti
! Date  : 28/Jul/2015
!

! Force (nxd,nzd) to be at most the product of a
! power of 2 and a single factor 3.
#define useFFTfit

MODULE dnsdata

  USE, intrinsic :: iso_c_binding
  USE rbmat
  USE mpi_transpose
  USE ffts
  IMPLICIT NONE

  !Simulation parameters
  real(C_DOUBLE) :: PI=3.1415926535897932384626433832795028841971
  integer(C_INT) :: nx,ny,nz,nxd,nzd
  real(C_DOUBLE) :: alfa0,beta0,ni,a,ymin,ymax,deltat,cflmax,time,dt_field,dt_save,t_max
  real(C_DOUBLE) :: meanpx,meanpz,meanflowx,meanflowz
  integer(C_INT), allocatable :: izd(:)
  complex(C_DOUBLE_COMPLEX), allocatable :: ialfa(:),ibeta(:)
  real(C_DOUBLE), allocatable :: k2(:,:)
  logical :: time_from_restart
  !Grid
  integer(C_INT), private :: iy
  real(C_DOUBLE), allocatable :: y(:),dy(:)
  real(C_DOUBLE) :: dx,dz,factor
  !Derivatives
  TYPE(Di), allocatable :: der(:)
  real(C_DOUBLE), dimension(-2:2) :: d040,d140,d14m1,d04n,d14n,d24n,d14np1
  real(C_DOUBLE), allocatable :: D0mat(:,:), etamat(:,:), D2vmat(:,:)
  !Fourier-transformable arrays (allocated in ffts.f90)
  complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:) :: VVdx, VVdz
  real(C_DOUBLE), pointer, dimension(:,:,:,:) :: rVVdx
  !Solution
  TYPE(RHSTYPE),  allocatable :: memrhs(:,:,:), oldrhs(:,:,:)
  complex(C_DOUBLE_COMPLEX), allocatable :: V(:,:,:,:)
  !Boundary conditions
  real(C_DOUBLE), dimension(-2:2) :: v0bc,v0m1bc,vnbc,vnp1bc,eta0bc,eta0m1bc,etanbc,etanp1bc
  TYPE(BCOND),    allocatable :: bc0(:,:), bcn(:,:)
  !Mean pressure correction
  real(C_DOUBLE), private :: corrpx=0.d0, corrpz=0.d0
  !ODE Library
  real(C_DOUBLE) :: RK1_rai(1:3)=(/ 120.0d0/32.0d0, 2.0d0, 0.0d0 /), &
                    RK2_rai(1:3)=(/ 120.0d0/8.0d0,  50.0d0/8.0d0,  34.0d0/8.0d0 /), &
                    RK3_rai(1:3)=(/ 120.0d0/20.0d0, 90.0d0/20.0d0, 50.0d0/20.0d0 /)
  !Outstats
  real(C_DOUBLE) :: cfl=0.0d0
  integer(C_SIZE_T) :: istep,nstep

  CONTAINS

  !--------------------------------------------------------------!
  !---------------------- Read input files ----------------------!
  SUBROUTINE read_dnsin()
    logical :: i
    OPEN(15, file='dns.in')
    READ(15, *) nx, ny, nz; READ(15, *) alfa0, beta0; nxd=3*(nx+1)/2;nzd=3*nz
#ifdef useFFTfit
    i=fftFIT(nxd); DO WHILE (.NOT. i); nxd=nxd+1; i=fftFIT(nxd); END DO
    i=fftFIT(nzd); DO WHILE (.NOT. i); nzd=nzd+1; i=fftFIT(nzd); END DO
#endif
    READ(15, *) ni; READ(15, *) a, ymin, ymax; ni=1/ni
    READ(15, *) meanpx, meanpz; READ(15, *) meanflowx, meanflowz
    READ(15, *) deltat, cflmax, time
    READ(15, *) dt_field, dt_save, t_max, time_from_restart
    READ(15, *) nstep
    CLOSE(15)
    dx=PI/(alfa0*nxd); dz=2.0d0*PI/(beta0*nzd);  factor=1.0d0/(2.0d0*nxd*nzd)
  END SUBROUTINE read_dnsin

  !--------------------------------------------------------------!
  !---------------- Allocate memory for solution ----------------!
  SUBROUTINE init_memory()
    INTEGER(C_INT) :: ix,iz
    ALLOCATE(V(-1:ny+1,-nz:nz,nx0:nxN,1:3))
    ALLOCATE(memrhs(0:2,-nz:nz,nx0:nxN),oldrhs(1:ny-1,-nz:nz,nx0:nxN),bc0(-nz:nz,nx0:nxN),bcn(-nz:nz,nx0:nxN))
#define newrhs(iy,iz,ix) memrhs(MOD(iy+1000,3),iz,ix)
#define imod(iy) MOD(iy+1000,5)
    ALLOCATE(der(1:ny-1),d0mat(1:ny-1,-2:2),etamat(1:ny-1,-2:2),D2vmat(1:ny-1,-2:2),y(-1:ny+1),dy(1:ny-1))
    ALLOCATE(izd(-nz:nz),ialfa(nx0:nxN),ibeta(-nz:nz),k2(-nz:nz,nx0:nxN))
    y=(/(ymin+0.5d0*(ymax-ymin)*(tanh(a*(2.0d0*real(iy)/real(ny)-1))/tanh(a)+0.5d0*(ymax-ymin)), iy=-1, ny+1)/)
    dy=(/( 0.5d0*(y(iy+1)-y(iy-1)) , iy=1, ny-1)/)
    izd=(/(merge(iz,nzd+iz,iz>=0),iz=-nz,nz)/);     ialfa=(/(dcmplx(0.0d0,ix*alfa0),ix=nx0,nxN)/);
    ibeta=(/(dcmplx(0.0d0,iz*beta0),iz=-nz,nz)/); 
    FORALL  (iz=-nz:nz,ix=nx0:nxN) k2(iz,ix)=(alfa0*ix)**2.0d0+(beta0*iz)**2.0d0
    IF (has_terminal) OPEN(UNIT=101,FILE='Runtimedata',ACTION='write')
  END SUBROUTINE init_memory

  !--------------------------------------------------------------!
  !--------------- Deallocate memory for solution ---------------!
  SUBROUTINE free_memory()
    DEALLOCATE(V,memrhs,oldrhs,der,bc0,bcn,d0mat,etamat,D2vmat,y,dy)
    IF (has_terminal) CLOSE(UNIT=101)
  END SUBROUTINE free_memory

  !--------------------------------------------------------------!
  !--------------- Set-up the compact derivatives ---------------!
  SUBROUTINE setup_derivatives()
    real(C_DOUBLE) :: M(0:4,0:4), t(0:4)
    integer(C_INT) :: iy,i,j
    DO iy=1,ny-1
      FORALL (i=0:4, j=0:4) M(i,j)=(y(iy-2+j)-y(iy))**(4.0d0-i); CALL LUdecomp(M)
      t=0; t(0)=24
      der(iy)%d4(-2:2)=M.bs.t
      FORALL (i=0:4, j=0:4) M(i,j)=(5.0d0-i)*(6.0d0-i)*(7.0d0-i)*(8.0d0-i)*(y(iy-2+j)-y(iy))**(4.0d0-i); CALL LUdecomp(M)
      FORALL (i=0:4) t(i)=sum( der(iy)%d4(-2:2)*(y(iy-2:iy+2)-y(iy))**(8.0d0-i) )
      der(iy)%d0(-2:2)=M.bs.t
      FORALL (i=0:4, j=0:4) M(i,j)=(y(iy-2+j)-y(iy))**(4.0d0-i); CALL LUdecomp(M)
      t=0; FORALL (i=0:2) t(i)=sum( der(iy)%d0(-2:2)*(4.0d0-i)*(3.0d0-i)*(y(iy-2:iy+2)-y(iy))**(2.0d0-i) )
      der(iy)%d2(-2:2)=M.bs.t
      t=0; FORALL (i=0:3) t(i)=sum( der(iy)%d0(-2:2)*(4.0d0-i)*(y(iy-2:iy+2)-y(iy))**(3.0d0-i) )
      der(iy)%d1(-2:2)=M.bs.t
    END DO
    FORALL (i=0:4, j=0:4) M(i,j)=(y(-1+j)-y(0))**(4.0d0-i); CALL LUdecomp(M)
    t=0; t(3)=1.0; d140(-2:2)=M.bs.t
    FORALL (i=0:4, j=0:4) M(i,j)=(y(-1+j)-y(-1))**(4.0d0-i); CALL LUdecomp(M)
    t=0; t(3)=1.0; d14m1(-2:2)=M.bs.t
    d04n=0; d04n(1)=1; d040=0; d040(-1)=1
    FORALL (i=0:4, j=0:4) M(i,j)=(y(ny-3+j)-y(ny))**(4.0d0-i); CALL LUdecomp(M)
    t=0; t(3)=1; d14n(-2:2)=M.bs.t
    t=0; t(2)=2; d24n(-2:2)=M.bs.t
    FORALL (i=0:4, j=0:4) M(i,j)=(y(ny-3+j)-y(ny+1))**(4.0d0-i); CALL LUdecomp(M)
    t=0; t(3)=1; d14np1(-2:2)=M.bs.t
    FORALL (iy=1:ny-1) D0mat(iy,-2:2)=der(iy)%d0(-2:2); CALL LU5decomp(D0mat)
  END SUBROUTINE setup_derivatives

  !--------------------------------------------------------------!
  !--------------- Set-up the boundary conditions ---------------!
  SUBROUTINE setup_boundary_conditions()
    v0bc=d040; v0m1bc=d140; eta0bc=d040
    vnbc=d04n; vnp1bc=d14n; etanbc=d04n
    etanp1bc=der(ny-1)%d4
    eta0m1bc=der(1)%d4
    v0bc(-1:2)=v0bc(-1:2)-v0bc(-2)*v0m1bc(-1:2)/v0m1bc(-2)
    eta0bc(-1:2)=eta0bc(-1:2)-eta0bc(-2)*eta0m1bc(-1:2)/eta0m1bc(-2)
    vnbc(-2:1)=vnbc(-2:1)-vnbc(2)*vnp1bc(-2:1)/vnp1bc(2)
    etanbc(-2:1)=etanbc(-2:1)-etanbc(2)*etanp1bc(-2:1)/etanp1bc(2)
  END SUBROUTINE setup_boundary_conditions

  !--------------------------------------------------------------!
  !---------------- integral in the y-direction -----------------!
  PURE FUNCTION yintegr(f) result(II)
    real(C_DOUBLE), intent(in) :: f(-1:ny+1)
    real(C_DOUBLE) :: II, yp1, ym1, a1, a2, a3
    integer(C_INT) :: iy
    II=0.0d0
    DO iy=1,ny-1,2
      yp1=y(iy+1)-y(iy); ym1=y(iy-1)-y(iy)
      a1=-1.0d0/3.0d0*ym1+1.0d0/6.0d0*yp1+1.0d0/6.0d0*yp1*yp1/ym1
      a3=+1.0d0/3.0d0*yp1-1.0d0/6.0d0*ym1-1.0d0/6.0d0*ym1*ym1/yp1
      a2=yp1-ym1-a1-a3
      II=II+a1*f(iy-1)+a2*f(iy)+a3*f(iy+1)
    END DO
  END FUNCTION yintegr

#define rD0(f,g,k) sum(dcmplx(der(iy)%d0(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d0(-2:2)*dreal(f(iy-2:iy+2,iz,ix,k))))
#define rD1(f,g,k) sum(dcmplx(der(iy)%d1(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d1(-2:2)*dreal(f(iy-2:iy+2,iz,ix,k))))
#define rD2(f,g,k) sum(dcmplx(der(iy)%d2(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d2(-2:2)*dreal(f(iy-2:iy+2,iz,ix,k))))
#define rD4(f,g,k) sum(dcmplx(der(iy)%d4(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d4(-2:2)*dreal(f(iy-2:iy+2,iz,ix,k))))
#define D0(f,g) sum(dcmplx(der(iy)%d0(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d0(-2:2)*dimag(f(iy-2:iy+2,iz,ix,g))))
#define D1(f,g) sum(dcmplx(der(iy)%d1(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d1(-2:2)*dimag(f(iy-2:iy+2,iz,ix,g))))
#define D2(f,g) sum(dcmplx(der(iy)%d2(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d2(-2:2)*dimag(f(iy-2:iy+2,iz,ix,g))))
#define D4(f,g) sum(dcmplx(der(iy)%d4(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d4(-2:2)*dimag(f(iy-2:iy+2,iz,ix,g))))
  !--------------------------------------------------------------!
  !---COMPLEX----- derivative in the y-direction ----------------!
  SUBROUTINE COMPLEXderiv(f0,f1)
    complex(C_DOUBLE_COMPLEX), intent(in)  :: f0(-1:ny+1)
    complex(C_DOUBLE_COMPLEX), intent(out) :: f1(-1:ny+1)
    f1(0)=sum(d140(-2:2)*f0(-1:3))
    f1(-1)=sum(d14m1(-2:2)*f0(-1:3))
    f1(ny)=sum(d14n(-2:2)*f0(ny-3:ny+1))
    f1(ny+1)=sum(d14np1(-2:2)*f0(ny-3:ny+1))
    DO CONCURRENT (iy=1:ny-1)
      f1(iy)=sum(der(iy)%d1(-2:2)*f0(iy-2:iy+2))
    END DO
    f1(1)=f1(1)-(der(1)%d0(-1)*f1(0)+der(1)%d0(-2)*f1(-1))
    f1(2)=f1(2)-der(2)%d0(-2)*f1(0)
    f1(ny-1)=f1(ny-1)-(der(ny-1)%d0(1)*f1(ny)+der(ny-1)%d0(2)*f1(ny+1))
    f1(ny-2)=f1(ny-2)-der(ny-2)%d0(2)*f1(ny)
    !f1(1:ny-1)=dcmplx(D0mat.bsr.dreal(f1(1:ny-1)),D0mat.bsr.dimag(f1(1:ny-1)))
    CALL LeftLU5div(D0mat,f1(1:ny-1))
  END SUBROUTINE COMPLEXderiv

  !--------------------------------------------------------------!
  !----------------- apply the boundary conditions --------------!
  PURE SUBROUTINE applybc_0(EQ,bc0,bc0m1)
    real(C_DOUBLE), intent(inout) :: EQ(1:ny-1,-2:2)
    real(C_DOUBLE), intent(in) :: bc0(-2:2),bc0m1(-2:2)
    EQ(1,-1:2)=EQ(1,-1:2)-EQ(1,-2)*bc0m1(-1:2)/bc0m1(-2)
    EQ(1, 0:2)=EQ(1, 0:2)-EQ(1,-1)*bc0(0:2)/bc0(-1)
    EQ(2,-1:1)=EQ(2,-1:1)-EQ(2,-2)*bc0(0:2)/bc0(-1)
  END SUBROUTINE applybc_0

  PURE SUBROUTINE applybc_n(EQ,bcn,bcnp1)
    real(C_DOUBLE), intent(inout) :: EQ(1:ny-1,-2:2)
    real(C_DOUBLE), intent(in) :: bcn(-2:2),bcnp1(-2:2)
    EQ(ny-1,-2:1)=EQ(ny-1,-2:1)-EQ(ny-1,2)*bcnp1(-2:1)/bcnp1(2)
    EQ(ny-1,-2:0)=EQ(ny-1,-2:0)-EQ(ny-1,1)*bcn(-2:0)/bcn(1)
    EQ(ny-2,-1:1)=EQ(ny-2,-1:1)-EQ(ny-2,2)*bcn(-2:0)/bcn(1)
  END SUBROUTINE applybc_n


#define OS(iy,j) (ni*(der(iy)%d4(j)-2.0d0*k2(iz,ix)*der(iy)%d2(j)+k2(iz,ix)*k2(iz,ix)*der(iy)%d0(j)))
#define SQ(iy,j) (ni*(der(iy)%d2(j)-k2(iz,ix)*der(iy)%d0(j)))
  !--------------------------------------------------------------!
  !------------------- solve the linear system  -----------------!
  SUBROUTINE linsolve(lambda)
    real(C_DOUBLE), intent(in) :: lambda
    integer(C_INT) :: ix,iz
    complex(C_DOUBLE_COMPLEX) :: temp(-1:ny+1)
    real(C_DOUBLE) :: ucor(-1:ny+1)
    DO iz=-nz,nz
      DO ix=nx0,nxN
        IF (ix==0 .AND. iz==0) THEN
          bc0(iz,ix)%v=0; bc0(iz,ix)%vy=0; bc0(iz,ix)%eta=dcmplx(dreal(bc0(iz,ix)%u)-dimag(bc0(iz,ix)%w),dimag(bc0(iz,ix)%u)+dreal(bc0(iz,ix)%w))
          bcn(iz,ix)%v=0; bcn(iz,ix)%vy=0; bcn(iz,ix)%eta=dcmplx(dreal(bcn(iz,ix)%u)-dimag(bcn(iz,ix)%w),dimag(bcn(iz,ix)%u)+dreal(bcn(iz,ix)%w))
        ELSE
          bc0(iz,ix)%vy=-ialfa(ix)*bc0(iz,ix)%u-ibeta(iz)*bc0(iz,ix)%w; bc0(iz,ix)%eta=ibeta(iz)*bc0(iz,ix)%u-ialfa(ix)*bc0(iz,ix)%w
          bcn(iz,ix)%vy=-ialfa(ix)*bcn(iz,ix)%u-ibeta(iz)*bcn(iz,ix)%w; bcn(iz,ix)%eta=ibeta(iz)*bcn(iz,ix)%u-ialfa(ix)*bcn(iz,ix)%w
        END IF
        bc0(iz,ix)%v=bc0(iz,ix)%v-v0bc(-2)*bc0(iz,ix)%vy/v0m1bc(-2)
        bcn(iz,ix)%v=bcn(iz,ix)%v-vnbc(2)*bcn(iz,ix)%vy/vnp1bc(2)
        DO CONCURRENT (iy=1:ny-1)
          D2vmat(iy,-2:2)=lambda*(der(iy)%d2(-2:2)-k2(iz,ix)*der(iy)%d0(-2:2))-OS(iy,-2:2)
          etamat(iy,-2:2)=lambda*der(iy)%d0(-2:2)-SQ(iy,-2:2)
        END DO
        CALL applybc_0(D2vmat,v0bc,v0m1bc)
        CALL applybc_n(D2vmat,vnbc,vnp1bc)
        V(1,iz,ix,2)=V(1,iz,ix,2)-D2vmat(1,-2)*bc0(iz,ix)%vy/v0m1bc(-2)-D2vmat(1,-1)*bc0(iz,ix)%v/v0bc(-1)
        V(2,iz,ix,2)=V(2,iz,ix,2)-D2vmat(2,-2)*bc0(iz,ix)%v/v0bc(-1)       
        V(ny-1,iz,ix,2)=V(ny-1,iz,ix,2)-D2vmat(ny-1,2)*bcn(iz,ix)%vy/vnp1bc(2)-D2vmat(ny-1,1)*bcn(iz,ix)%v/vnbc(1)
        V(ny-2,iz,ix,2)=V(ny-2,iz,ix,2)-D2vmat(ny-2,2)*bcn(iz,ix)%v/vnbc(1)
        CALL applybc_0(etamat,eta0bc,eta0m1bc)
        CALL applybc_n(etamat,etanbc,etanp1bc)
        V(1,iz,ix,1)=V(1,iz,ix,1)-etamat(1,-1)*bc0(iz,ix)%eta/eta0bc(-1)
        V(2,iz,ix,1)=V(2,iz,ix,1)-etamat(2,-2)*bc0(iz,ix)%eta/eta0bc(-1)
        V(ny-1,iz,ix,1)=V(ny-1,iz,ix,1)-etamat(ny-1,1)*bcn(iz,ix)%eta/etanbc(1)
        V(ny-2,iz,ix,1)=V(ny-2,iz,ix,1)-etamat(ny-2,2)*bcn(iz,ix)%eta/etanbc(1)
        CALL LU5decomp(D2vmat); CALL LU5decomp(etamat)
        CALL LeftLU5div(D2vmat,V(1:ny-1,iz,ix,2))
        V(0,iz,ix,2)=(bc0(iz,ix)%v-sum(V(1:3,iz,ix,2)*v0bc(0:2)))/v0bc(-1)
        V(-1,iz,ix,2)=(bc0(iz,ix)%vy-sum(V(0:3,iz,ix,2)*v0m1bc(-1:2)))/v0m1bc(-2)
        V(ny,iz,ix,2)=(bcn(iz,ix)%v-sum(V(ny-3:ny-1,iz,ix,2)*vnbc(-2:0)))/vnbc(1)
        V(ny+1,iz,ix,2)=(bcn(iz,ix)%vy-sum(V(ny-3:ny,iz,ix,2)*vnp1bc(-2:1)))/vnp1bc(2)
        CALL LeftLU5div(etamat,V(1:ny-1,iz,ix,1))
        V(0,iz,ix,1)=(bc0(iz,ix)%eta-sum(V(1:3,iz,ix,1)*eta0bc(0:2)))/eta0bc(-1)
        V(-1,iz,ix,1)=-sum(V(0:3,iz,ix,1)*eta0m1bc(-1:2))/eta0m1bc(-2)
        V(ny,iz,ix,1)=(bcn(iz,ix)%eta-sum(V(ny-3:ny-1,iz,ix,1)*etanbc(-2:0)))/etanbc(1)
        V(ny+1,iz,ix,1)=-sum(V(ny-3:ny,iz,ix,1)*etanp1bc(-2:1))/etanp1bc(2)
        IF (ix==0 .AND. iz==0) THEN
            V(:,0,0,3) = dcmplx(dimag(V(:,0,0,1)),0.d0); 
            V(:,0,0,1) = dcmplx(dreal(V(:,0,0,1)),0.d0); 
            ucor(-1:0)=0; ucor(1:ny-1)=1; ucor(ny:ny+1)=0
            ucor(1:ny-1)=etamat.bsr.ucor(1:ny-1)
            ucor(0)=-sum(ucor(1:3)*eta0bc(0:2))/eta0bc(-1)
            ucor(-1)=-sum(ucor(0:3)*eta0m1bc(-1:2))/eta0m1bc(-2)
            ucor(ny)=-sum(ucor(ny-3:ny-1)*etanbc(-2:0))/etanbc(1)
            ucor(ny+1)=-sum(ucor(ny-3:ny)*etanp1bc(-2:1))/etanp1bc(2)          
            IF (abs(meanflowx)>1.0d-7) THEN
              corrpx = (meanflowx-yintegr(dreal(V(:,0,0,1))))/yintegr(ucor)
              V(:,0,0,1)=dcmplx(dreal(V(:,0,0,1))+corrpx*ucor,dimag(V(:,0,0,1)))
            END IF
            IF (abs(meanflowz)>1.0d-7) THEN
              corrpz = (meanflowz-yintegr(dreal(V(:,0,0,3))))/yintegr(ucor)
              V(:,0,0,3)=dcmplx(dreal(V(:,0,0,3))+corrpz*ucor,dimag(V(:,0,0,3)))
            END IF
        ELSE
            CALL COMPLEXderiv(V(:,iz,ix,2),V(:,iz,ix,3))
            temp=(ialfa(ix)*V(:,iz,ix,3)-ibeta(iz)*V(:,iz,ix,1))/k2(iz,ix)
            V(:,iz,ix,3)=(ibeta(iz)*V(:,iz,ix,3)+ialfa(ix)*V(:,iz,ix,1))/k2(iz,ix)
            V(:,iz,ix,1)=temp
        END IF
      END DO
    END DO
  END SUBROUTINE linsolve

   !--------------------------------------------------------------!
   !------------------------ convolutions ------------------------!
  SUBROUTINE convolutions(iy,i,compute_cfl)
     integer(C_INT), intent(in) :: iy,i
     logical, intent(in) :: compute_cfl
     integer(C_INT) :: ix,iz,iV
     VVdz(1:nz+1,1:nxB,1:3,i)=V(iy,0:nz,nx0:nxN,1:3);         VVdz(nz+2:nzd-nz,1:nxB,1:3,i)=0;
     VVdz(nzd+1-nz:nzd,1:nxB,1:3,i)=V(iy,-nz:-1,nx0:nxN,1:3); 
     DO iV=1,3
       CALL IFT(VVdz(1:nzd,1:nxB,iV,i)); CALL MPI_Alltoall(VVdz(:,:,iV,i), 1, Mdz, VVdx(:,:,iV,i), 1, Mdx, MPI_COMM_WORLD)
       VVdx(nx+2:nxd+1,1:nzB,iV,i)=0;    CALL RFT(VVdx(1:nxd+1,1:nzB,iV,i),rVVdx(1:2*nxd+2,1:nzB,iV,i));
     END DO
     VVdx(nx+2:nxd+1,1:nzB,4:6,i)=0;
     IF (compute_cfl .and. iy>=1 .and. iy<=ny-1) THEN
           cfl=max(cfl,(maxval(abs(rVVdx(1:2*nxd,1:nzB,1,i))/dx     + &
                               abs(rVVdx(1:2*nxd,1:nzB,2,i))/dy(iy) + &
                               abs(rVVdx(1:2*nxd,1:nzB,3,i))/dz)))
     END IF
     rVVdx(1:2*nxd,1:nzB,4,i)  = rVVdx(1:2*nxd,1:nzB,1,i)  * rVVdx(1:2*nxd,1:nzB,2,i)*factor
     rVVdx(1:2*nxd,1:nzB,5,i)  = rVVdx(1:2*nxd,1:nzB,2,i)  * rVVdx(1:2*nxd,1:nzB,3,i)*factor
     rVVdx(1:2*nxd,1:nzB,6,i)  = rVVdx(1:2*nxd,1:nzB,1,i)  * rVVdx(1:2*nxd,1:nzB,3,i)*factor
     rVVdx(1:2*nxd,1:nzB,1:3,i)= rVVdx(1:2*nxd,1:nzB,1:3,i)* rVVdx(1:2*nxd,1:nzB,1:3,i)*factor
     DO iV=1,6
       CALL HFT(rVVdx(1:2*nxd+2,1:nzB,iV,i),VVdx(1:nxd+1,1:nzB,iV,i)); 
       CALL MPI_Alltoall(VVdx(:,:,iV,i), 1, Mdx, VVdz(:,:,iV,i), 1, Mdz, MPI_COMM_WORLD)
       CALL FFT(VVdz(1:nzd,1:nxB,iV,i));
     END DO
   END SUBROUTINE convolutions


  !--------------------------------------------------------------!
  !-------------------------- buildRHS --------------------------!
  ! (u,v,w) = (1,2,3)
  ! (uu,vv,ww,uv,vw,uw) = (1,2,3,4,5,6)
#define DD(f,k) ( der(iy)%f(-2)*VVdz(izd(iz)+1,ix+1-nx0,k,im2)+der(iy)%f(-1)*VVdz(izd(iz)+1,ix+1-nx0,k,im1)+der(iy)%f(0)*VVdz(izd(iz)+1,ix+1-nx0,k,i0)+ \
                  der(iy)%f(1 )*VVdz(izd(iz)+1,ix+1-nx0,k,i1 )+der(iy)%f(2 )*VVdz(izd(iz)+1,ix+1-nx0,k,i2 ) )
#define timescheme(rhs,old,unkn,impl,expl) rhs=ODE(1)*(unkn)/deltat+(impl)+ODE(2)*(expl)-ODE(3)*(old); old=expl                   
  SUBROUTINE buildrhs(ODE,compute_cfl)
    logical, intent(in) :: compute_cfl
    real(C_DOUBLE), intent(in) :: ODE(1:3)
    integer(C_INT) :: iy,iz,ix,im2,im1,i0,i1,i2
    complex(C_DOUBLE_COMPLEX) :: rhsu,rhsv,rhsw,DD0_6,DD1_6,expl
    DO iy=-3,ny+1
      IF (iy<=ny-1) THEN
      CALL convolutions(iy+2,imod(iy+2)+1,compute_cfl)
      IF (iy>=1) THEN
        im2=imod(iy-2)+1; im1=imod(iy-1)+1; i0=imod(iy)+1; i1=imod(iy+1)+1; i2=imod(iy+2)+1;
        DO iz=-nz,nz 
        DO ix=nx0,nxN
            DD0_6=DD(d0,6); DD1_6=DD(d1,6);
            rhsu=-ialfa(ix)*DD(d0,1)-DD(d1,4)-ibeta(iz)*DD0_6
            rhsv=-ialfa(ix)*DD(d0,4)-DD(d1,2)-ibeta(iz)*DD(d0,5)
            rhsw=-ialfa(ix)*DD0_6-DD(d1,5)-ibeta(iz)*DD(d0,3)
            expl=(ialfa(ix)*(ialfa(ix)*DD(d1,1)+DD(d2,4)+ibeta(iz)*DD1_6)+&
                  ibeta(iz)*(ialfa(ix)*DD1_6+DD(d2,5)+ibeta(iz)*DD(d1,3))-k2(iz,ix)*rhsv &
                 )
            timescheme(newrhs(iy,iz,ix)%D2v, oldrhs(iy,iz,ix)%D2v, D2(V,2)-k2(iz,ix)*D0(V,2),sum(OS(iy,-2:2)*V(iy-2:iy+2,iz,ix,2)),expl); !(D2v)
            IF (ix==0 .AND. iz==0) THEN
              expl=(dcmplx(dreal(rhsu)+meanpx,dreal(rhsw)+meanpz) &
                   )
              timescheme(newrhs(iy,0,0)%eta,oldrhs(iy,0,0)%eta,rD0(V,1,3),ni*rD2(V,1,3),expl) !(Ubar, Wbar)
            ELSE
              expl=(ibeta(iz)*rhsu-ialfa(ix)*rhsw &
                   )
              timescheme(newrhs(iy,iz,ix)%eta, oldrhs(iy,iz,ix)%eta,ibeta(iz)*D0(V,1)-ialfa(ix)*D0(V,3),sum(SQ(iy,-2:2)*[ibeta(iz)*V(iy-2:iy+2,iz,ix,1)-ialfa(ix)*V(iy-2:iy+2,iz,ix,3)]),expl) !(eta)
            END IF
        END DO
        END DO
      END IF
      END IF
      IF (iy-2>=1) THEN
        DO CONCURRENT (ix=nx0:nxN, iz=-nz:nz) 
          V(iy-2,iz,ix,1) = newrhs(iy-2,iz,ix)%eta; V(iy-2,iz,ix,2) = newrhs(iy-2,iz,ix)%d2v; 
        END DO
      END IF      
    END DO
  END SUBROUTINE buildrhs

  !--------------------------------------------------------------!
  !-------------------- read_restart_file -----------------------! 
  SUBROUTINE read_restart_file()
    integer(C_SIZE_T) :: iV,ix,iy,iz,io,nxB_t,nx_t,nz_t,ny_t,iproc_t,br=8,bc=16,iV_t,b1=1,b7=7,b3=3
    integer(C_SIZE_T) :: pos
    OPEN(UNIT=100,FILE="Dati.cart.out",access="stream",status="old",action="read",iostat=io)
    nx_t=nx+1; ny_t=ny+3; nz_t=2*nz+1; iproc_t=iproc; nxB_t=nxB
    IF (io==0) THEN
      IF (has_terminal) WRITE(*,*) "Reading restart file..."
      READ(100,POS=1) nx,ny,nz,alfa0,beta0,ni,a,ymin,ymax,time
      DO iV=1,3
        pos=bc*ny_t*nz_t*nxB_t*iproc_t+(iV-b1)*(bc*ny_t*nz_t*nx_t)+b1+(br*b7+b3*SIZEOF(nx))
        WRITE(*,*) pos,iproc
        READ(100,POS=pos) V(:,:,:,iV)
      END DO
      CLOSE(100)
    ELSE
      V=0
      IF (has_terminal) WRITE(*,*) "Generating initial field..."
      DO iy=-1,ny+1; DO ix=nx0,nxN; DO iz=-nz,nz
          V(iy,iz,ix,1) = 0.0001*EXP(dcmplx(0,RAND()-0.5));  V(iy,iz,ix,2) = 0.0001*EXP(dcmplx(0,RAND()-0.5));  V(iy,iz,ix,3) = 0.0001*EXP(dcmplx(0,RAND()-0.5));
      END DO;        END DO;        END DO
      IF (has_terminal) THEN
        DO CONCURRENT (iy=-1:ny+1)
          V(iy,0,0,1)=y(iy)*(2-y(iy))*3.d0/2.d0 + 0.001*SIN(8*y(iy)*2*PI);
        END DO
      END IF
    END IF
  END SUBROUTINE read_restart_file

  !--------------------------------------------------------------!
  !-------------------- save_restart_file -----------------------!
  SUBROUTINE save_restart_file(filename)
    integer(C_SIZE_T) :: iV,ix,iy,iz,io,nxB_t,nx_t,nz_t,ny_t,iproc_t,br=8,bc=16,iV_t,b1=1,b7=7,b3=3
    integer(C_SIZE_T) :: pos,i
    character(len=40) :: filename
    DO i=0,nproc-1
      IF (i==iproc) THEN
        OPEN(UNIT=100,FILE=TRIM(filename),access="stream",action="write")
        nx_t=nx+1; ny_t=ny+3; nz_t=2*nz+1; iproc_t=iproc; nxB_t=nxB
        IF (has_terminal) WRITE(UNIT=100,POS=1) nx,ny,nz,alfa0,beta0,ni,a,ymin,ymax,time
        DO iV=1,3
          pos=bc*ny_t*nz_t*nxB_t*iproc_t+(iV-b1)*(bc*ny_t*nz_t*nx_t)+b1+(br*b7+b3*SIZEOF(nx))
          WRITE(100,POS=pos) V(:,:,:,iV)
        END DO
        CLOSE(100)
      END IF
      CALL MPI_Barrier(MPI_COMM_WORLD)
    END DO
  END SUBROUTINE save_restart_file


  !--------------------------------------------------------------!
  !------------------------- outstats ---------------------------!
  SUBROUTINE outstats()
   real(C_DOUBLE) :: runtime_global   !cfl
   character(len=40) :: istring, filename
   CALL MPI_Allreduce(cfl,runtime_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD); cfl=0
   IF (cflmax>0)  deltat=cflmax/runtime_global;
   IF (has_terminal) THEN
     WRITE(*,"(F6.4,3X,4(F11.6,3X),4(F9.4,3X),2(F9.6,3X))") &
           time,sum(d140(-2:2)*dreal(V(-1:3,0,0,1))),-sum(d14n(-2:2)*dreal(V(ny-3:ny+1,0,0,1))),&
                sum(d140(-2:2)*dreal(V(-1:3,0,0,3))),-sum(d14n(-2:2)*dreal(V(ny-3:ny+1,0,0,3))),&
                yintegr(dreal(V(:,0,0,1))),meanpx+corrpx,yintegr(dreal(V(:,0,0,3))),meanpz +corrpz,&
                runtime_global*deltat,deltat
     WRITE(101,*) time,sum(d140(-2:2)*dreal(V(-1:3,0,0,1))),-sum(d14n(-2:2)*dreal(V(ny-3:ny+1,0,0,1))),&
                       sum(d140(-2:2)*dreal(V(-1:3,0,0,3))),-sum(d14n(-2:2)*dreal(V(ny-3:ny+1,0,0,3))),&
                       yintegr(dreal(V(:,0,0,1))),meanpx+corrpx,yintegr(dreal(V(:,0,0,3))),meanpz +corrpz,&
                       runtime_global*deltat,deltat
   END IF
   runtime_global=0
   !Save Dati.cart.out
   IF ( ((FLOOR((time+0.5*deltat)/dt_save) > FLOOR((time-0.5*deltat)/dt_save)) .AND. (istep>1)) .OR. istep==nstep ) THEN
     IF (has_terminal) WRITE(*,*) "Writing Dati.cart.out at time ", time
     filename="Dati.cart.out"; CALL save_restart_file(filename)
   END IF
   IF ( (FLOOR((time+0.5*deltat)/dt_field) > FLOOR((time-0.5*deltat)/dt_field)) .AND. (time>0) ) THEN
     WRITE(istring,*) FLOOR(time/dt_field)
     IF (has_terminal) WRITE(*,*) "Writing Dati.cart."//TRIM(ADJUSTL(istring))//".out at time ", time
     filename="Dati.cart."//TRIM(ADJUSTL(istring))//".out"; CALL save_restart_file(filename)
   END IF
  END SUBROUTINE outstats

END MODULE dnsdata
