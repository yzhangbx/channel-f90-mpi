!============================================!
!                                            !
!             Distributed Traspose           !
!                  for the                   !
!      Direct Numerical Simulation (DNS)     !
!        of a turbulent channel flow         !
!                                            !
!============================================!
! 
! Author: Dr. Davide Gatti
! Date  : 23/Sep/2015
! 

MODULE mpi_transpose
  
  USE, intrinsic :: iso_c_binding
  USE typedef
  USE mpi_f08
  IMPLICIT NONE

  integer(C_INT),save :: nproc,iproc,ierr
  integer(C_INT), save :: nx0,nxN,nxB,nz0,nzN,nzB,block
  TYPE(MPI_Datatype), save :: Mdz3,Mdx3,Mdx6,Mdz6,cmpl,vel,momfl
  logical, save :: has_terminal
  

CONTAINS

  !------- Divide the problem in 1D slices -------! 
  !-----------------------------------------------!
  SUBROUTINE init_MPI(iproc,nproc,nx,nxd,nzd,nx0,nxN,nxB,nz0,nzN,nzB,block)
    integer(C_INT), intent(in)  :: nx,nxd,nzd,iproc,nproc
    integer(C_INT), intent(out) :: nx0,nxN,nxB,nz0,nzN,nzB,block
    TYPE(MPI_Datatype) :: row,column,tmp
    integer(kind=MPI_ADDRESS_KIND) :: stride,lb
    ! Define which process write on screen
    has_terminal=(iproc==0)
    ! Calculate domain division (XXX add checks XXX)
    nx0=iproc*nx/nproc;   nxN=(iproc+1)*nx/nproc-1;  nxB=nxN-nx0+1;
    nz0=iproc*nzd/nproc;  nzN=(iproc+1)*nzd/nproc-1; nzB=nzN-nz0+1;
    block=max(nxB*nzd,nx*nzB)
    WRITE(*,*) "iproc=",iproc," nx0=",nx0," nxN=",nxN," nxB=",nxB, "nz0=",nz0," nzN=",nzN," nzB=",nzB, " s=", block 
    ! MPI derived datatyped - basics
    CALL MPI_Type_contiguous(2,MPI_DOUBLE_PRECISION,cmpl)   !complex
    CALL MPI_Type_contiguous(3,cmpl,vel); tmp=vel           !velocity
    lb=0; stride=6*16; CALL MPI_Type_create_resized(tmp,lb,stride,vel);
    CALL MPI_Type_contiguous(6,cmpl,momfl)                  !momflux
    CALL MPI_Type_commit(cmpl); CALL MPI_Type_commit(vel);
    CALL MPI_Type_commit(momfl);
    ! interlaved MPI datatypes - 3 variables
    CALL MPI_Type_vector(nxB,nzB,nzd,vel,row)
    lb=0; stride=6*16*nzB; CALL MPI_Type_create_resized(row,lb,stride,Mdz3);
    CALL MPI_Type_commit(Mdz3)    
    CALL MPI_Type_vector(nzB,1,nxd,vel,column)
    lb=0; stride=6*16;  CALL MPI_Type_create_resized(column,lb,stride,tmp)
    CALL MPI_Type_contiguous(nxB,tmp,Mdx3);
    CALL MPI_Type_commit(Mdx3)    
    ! interlaved MPI datatypes - 6 variables
    CALL MPI_Type_vector(nxB,nzB,nzd,momfl,row)
    lb=0; stride=6*16*nzB; CALL MPI_Type_create_resized(row,lb,stride,Mdz6);
    CALL MPI_Type_commit(Mdz6)    
    CALL MPI_Type_vector(nzB,1,nxd,momfl,column)
    lb=0; stride=6*16;  CALL MPI_Type_create_resized(column,lb,stride,tmp)
    CALL MPI_Type_contiguous(nxB,tmp,Mdx6);
    CALL MPI_Type_commit(Mdx6)    
  END SUBROUTINE init_MPI
  

END MODULE mpi_transpose
