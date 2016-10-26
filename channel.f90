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
#define chron


PROGRAM channel

  USE dnsdata
#ifdef crhon
  REAL timei,timee 
#endif

  ! Init MPI
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  CALL read_dnsin()
  CALL init_MPI(iproc,nproc,nx+1,nz,ny,nxd+1,nzd,nx0,nxN,nxB,nz0,nzN,nzB,block)
  CALL init_memory()

  ! Init various subroutines
  CALL init_fft(VVdz,VVdx,rVVdx,nxd,nx0,nxN,nxB,nzd,nz0,nzN,nzB)
  CALL setup_derivatives()
  CALL setup_boundary_conditions()
  CALL read_restart_file()

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
  END DO timeloop

  ! Realease memory
   CALL free_fft()
   CALL free_memory()
   CALL MPI_Finalize()


END PROGRAM channel
