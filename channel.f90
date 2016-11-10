!============================================!
!                                            !
!     Direct Numerical Simulation (DNS)      !
!       of a turbulent channel flow          !
!                                            !
!============================================!
!
! This program has been written following the
! KISS (Keep it Simple and Stupid) philosophy
!
! Author: Dr.-Ing. Davide Gatti
! Date  : 28/Jul/2015
!

! Measure per timestep execution time
!#define chron

PROGRAM channel

  USE dnsdata
#ifdef crhon
  REAL timei,timee
#endif
  ! Stats
  INTEGER(C_SIZE_T) :: istep, nstats=20, istats=0
  LOGICAL :: io
  REAL(C_DOUBLE), allocatable :: localstats(:,:), globalstats(:,:)
  REAL(C_DOUBLE) :: c

  ! Init MPI
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  CALL read_dnsin()
  CALL init_MPI(iproc,nproc,nx+1,nz,ny,nxd+1,nzd,nx0,nxN,nxB,nz0,nzN,nzB,block)
  CALL init_memory()

  ! Init various subroutines
  CALL init_fft(VVdz,VVdx,rVVdx,nxd,nxB,nzd,nzB)
  CALL setup_derivatives()
  CALL setup_boundary_conditions()
  CALL read_restart_file()

  ! Allocate memory for stats
  ALLOCATE(localstats(-1:ny+1,1:5),globalstats(-1:ny+1,1:5))
  IF (has_terminal) THEN 
    INQUIRE(FILE="stats.dat",EXIST=io)
    IF (io==0) THEN 
      istats=0; globalstats=0
      OPEN(UNIT=102,FILE='stats.dat',STATUS="new",ACCESS='stream',ACTION='readwrite')
    ELSE
      OPEN(UNIT=102,FILE='stats.dat',STATUS="old",ACCESS='stream',ACTION='readwrite')
      READ(102,POS=1) istats, globalstats
      globalstats=globalstats*istats
    END IF
  END IF

IF (has_terminal) THEN
  ! Output DNS.in
  WRITE(*,*) " "
  WRITE(*,*) "!====================================================!"
  WRITE(*,*) "!                     D   N   S                      !"
  WRITE(*,*) "!====================================================!"
  WRITE(*,*) " "
  WRITE(*,"(A,I5,A,I5,A,I5)") "   nx =",nx,"   ny =",ny,"   nz =",nz
  WRITE(*,"(A,I5,A,I5)") "   nxd =",nxd,"  nzd =",nzd
  WRITE(*,"(A,F6.4,A,F6.4,A,F8.6)") "   alfa0 =",alfa0,"       beta0 =",beta0,"   ni =",ni
  WRITE(*,"(A,F6.4,A,F6.4)") "   meanpx =",meanpx,"      meanpz =",meanpz
  WRITE(*,"(A,F6.4,A,F6.4)") "   meanflowx =",meanflowx, "   meanflowz =", meanflowz
  WRITE(*,*) " "
END IF

  ! Compute CFL
  DO iy=1,ny-1
   CALL convolutions(iy,1,.TRUE.)
  END DO
  ! Time loop
  CALL outstats()
  timeloop: DO WHILE (time<t_max-deltat/2.0)
#ifdef chron
    CALL CPU_TIME(timei)
#endif
    time=time+2.0/RK1_rai(1)*deltat
    CALL buildrhs(RK1_rai,.FALSE. ); CALL linsolve(RK1_rai(1)/deltat)
    time=time+2.0/RK2_rai(1)*deltat
    CALL buildrhs(RK2_rai,.FALSE.); CALL linsolve(RK2_rai(1)/deltat)
    time=time+2.0/RK3_rai(1)*deltat
    CALL buildrhs(RK3_rai,.TRUE.); CALL linsolve(RK3_rai(1)/deltat)
    CALL outstats()
#ifdef chron
    CALL CPU_TIME(timee)
    IF (has_terminal) WRITE(*,*) timee-timei
#endif
    ! Compute statistics
    istep=istep+1
    IF (MOD(istep,nstats)==0) THEN
      localstats=0; istats=istats+1
      IF (has_terminal) localstats(:,1)=localstats(:,1)+dreal(V(:,0,0,1))  ! U
      DO ix=nx0,nxN
        c = MERGE(1.0d0,2.0d0,ix==0) 
        localstats(:,5) = localstats(:,5) +c*SUM(V(:,:,ix,1)*CONJG(V(:,:,ix,2)),2)  ! uv
        localstats(:,2) = localstats(:,2) +c*SUM(V(:,:,ix,1)*CONJG(V(:,:,ix,1)),2)  ! uu
        localstats(:,3) = localstats(:,3) +c*SUM(V(:,:,ix,2)*CONJG(V(:,:,ix,2)),2)  ! vv
        localstats(:,4) = localstats(:,4) +c*SUM(V(:,:,ix,3)*CONJG(V(:,:,ix,3)),2)  ! ww
      END DO
      IF (has_terminal) THEN
        CALL MPI_Reduce(MPI_IN_PLACE,localstats,(ny+3)*5,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD)
      ELSE 
        CALL MPI_Reduce(localstats,localstats,(ny+3)*5,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD)
      END IF
      IF (has_terminal) THEN
        globalstats=globalstats+localstats; WRITE(102,POS=1) istats,globalstats/istats
      END IF
    END IF
  END DO timeloop

  IF (has_terminal) CLOSE(102)
  ! Realease memory
   CALL free_fft()
   CALL free_memory()
   CALL MPI_Finalize()


END PROGRAM channel
