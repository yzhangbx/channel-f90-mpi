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
  TYPE(MPI_Datatype), save :: Mdz,Mdx,cmpl,iVc,vel,velc
  logical, save :: has_terminal

CONTAINS

  !------- Divide the problem in 1D slices -------! 
  !-----------------------------------------------!
  SUBROUTINE init_MPI(iproc,nproc,nx,nz,ny,nxd,nzd,nx0,nxN,nxB,nz0,nzN,nzB,block)
    integer(C_INT), intent(in)  :: nx,nz,ny,nxd,nzd,iproc,nproc
    integer(C_INT), intent(out) :: nx0,nxN,nxB,nz0,nzN,nzB,block
    TYPE(MPI_Datatype) :: row,column,tmp
    integer(kind=MPI_ADDRESS_KIND) :: stride,lb
    ! Define which process write on screen
    has_terminal=(iproc==0)
    ! Calculate domain division (XXX add checks XXX)
    nx0=iproc*(nx)/nproc; nxN=(iproc+1)*(nx)/nproc-1;  nxB=nxN-nx0+1;
    nz0=iproc*nzd/nproc;  nzN=(iproc+1)*nzd/nproc-1; nzB=nzN-nz0+1;
    block=max(nxB*nzd,nx*nzB)
    WRITE(*,*) "iproc=",iproc," nx0=",nx0," nxN=",nxN," nxB=",nxB, "nz0=",nz0," nzN=",nzN," nzB=",nzB, " s=", block 
    ! MPI derived datatyped - basics
    CALL MPI_Type_contiguous(2,MPI_DOUBLE_PRECISION,cmpl)   !complex
    CALL MPI_Type_commit(cmpl)
    ! interlaved MPI datatypes - communicate plane of data
    CALL MPI_Type_vector(nxB,nzB,nzd,cmpl,row)
    lb=0; stride=8*2*nzB; CALL MPI_Type_create_resized(row,lb,stride,Mdz)
    CALL MPI_Type_commit(Mdz)    
    CALL MPI_Type_vector(nzB,1,nxd,cmpl,column)
    lb=0; stride=8*2;  CALL MPI_Type_create_resized(column,lb,stride,tmp)
    CALL MPI_Type_contiguous(nxB,tmp,Mdx)
    CALL MPI_Type_commit(Mdx)    
    ! interlaved MPI datatypes - map velocity field to file
    CALL MPI_Type_contiguous(nxB*(2*nz+1)*(ny+3),cmpl,iVc);
    CALL MPI_Type_commit(iVc)
    CALL MPI_Type_vector(3,1,nproc,iVc,vel); tmp=vel
    lb=8*2*nxB*(2*nz+1)*(ny+3)*iproc; stride=8*2*nx*(2*nz+1)*(ny+3); 
    CALL MPI_Type_create_resized(tmp,lb,stride,vel)
    CALL MPI_Type_commit(vel)
  END SUBROUTINE init_MPI
  

END MODULE mpi_transpose
