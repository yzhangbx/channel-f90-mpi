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
! Author: M.Sc. Davide Gatti
! Date  : 28/Jul/2015
! 

! MISSING: free-to-choose boundary conditions
!          computation of the CFL and other statistics

PROGRAM channel

  USE mpi_f08
  USE dnsdata
 
  ! Init MPI
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  CALL init_MPI(iproc,nproc,nx+1,nxd+1,nzd,nx0,nxN,nxB,nz0,nzN,nzB,block)
  CALL init_memory()

  ! Init various subroutines
  CALL init_fft(VVdz,VVdx,rVVdx,nxd,nxB,nzd,nzB)
  CALL setup_derivatives()
  CALL setup_boundary_conditions()
  CALL read_restart_file(V)

  ! Output DNS.in
  IF (has_terminal) THEN
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

  DO iy=2,ny-1
   CALL convolutions(iy,0,.TRUE.)
  END DO
  !Time loop
  CALL outstats()
  timeloop: DO WHILE (time<t_max-deltat/2.0) 
    time=time+2.0/RK1_rai_coeff*deltat
    CALL buildrhs(RK1_rai,.TRUE. ); CALL linsolve(RK1_rai_coeff/deltat)
    time=time+2.0/RK2_rai_coeff*deltat
    CALL buildrhs(RK2_rai,.FALSE.); CALL linsolve(RK2_rai_coeff/deltat)
    time=time+2.0/RK3_rai_coeff*deltat
    CALL buildrhs(RK3_rai,.FALSE.); CALL linsolve(RK3_rai_coeff/deltat)
    CALL outstats()
   END DO timeloop

  !Realease memory
  CALL free_fft()
  CALL free_memory()

END PROGRAM channel
