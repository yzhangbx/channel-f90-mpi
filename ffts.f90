!============================================!
!                                            !
!           Fast Fourier Transforms          !
!                  for the                   !
!      Direct Numerical Simulation (DNS)     !
!        of a turbulent channel flow         !
!                                            !
!============================================!
! 
! Author: M.Sc. Davide Gatti
! Date  : 28/Jul/2015
! 

MODULE ffts

  USE, intrinsic :: iso_c_binding
  USE typedef
  IMPLICIT NONE
  INCLUDE 'fftw3.f03'

  integer, save        :: plan_type=FFTW_PATIENT
  TYPE(C_PTR), save    :: pFFT,pIFT,pRFT,pHFT,ptrVVdx,ptrVVdz

CONTAINS

  !--------------------------------------------------------------!
  !------------------ init aligned FFT vectors ------------------!
  SUBROUTINE init_fft(VVdz,VVdx,rVVdx,nxd,nxB,nzd,nzB)
    integer(C_INT), intent(in) :: nxd,nxB,nzd,nzB
    complex(C_DOUBLE_COMPLEX), pointer, intent(out) :: VVdz(:,:,:)
    complex(C_DOUBLE_COMPLEX), pointer, intent(out) :: VVdx(:,:,:)
    real(C_DOUBLE), pointer, intent(out) :: rVVdx(:,:,:)
    integer(C_INT), dimension(1) :: n_z, n_x, rn_x
    n_z=[nzd]; n_x=[nxd]; rn_x=[2*nxd];
    ptrVVdz=fftw_alloc_complex(int(nxB*nzd*6, C_SIZE_T))
    ptrVVdx=fftw_alloc_complex(int((nxd+1)*nzB*6, C_SIZE_T))
    CALL c_f_pointer(ptrVVdz, VVdz, [6,nzd,nxB]);  CALL c_f_pointer(ptrVVdx,  VVdx, [6,nxd+1,nzB])
                                                   CALL c_f_pointer(ptrVVdx, rVVdx, [6,2*(nxd+1),nzB])    
    pFFT=fftw_plan_many_dft(1, n_z, 6, VVdz(1:6,:,1), n_z, 6, 1, VVdz(1:6,:,1), n_z, 6, 1, FFTW_FORWARD,  plan_type)
    pIFT=fftw_plan_many_dft(1, n_z, 3, VVdz(1:6,:,1), n_z, 6, 1, VVdz(1:6,:,1), n_z, 6, 1, FFTW_BACKWARD, plan_type)
    pRFT=fftw_plan_many_dft_c2r(1, rn_x, 3, VVdx(1:6,:,1),  n_x+1, 6, 1, rVVdx(1:6,:,1), rn_x+2,  6, 1,   plan_type)
    pHFT=fftw_plan_many_dft_r2c(1, rn_x, 6, rVVdx(1:6,:,1), rn_x+2, 6, 1, VVdx(1:6,:,1),  n_x+1,  6, 1,   plan_type)
  END SUBROUTINE init_fft

  SUBROUTINE FFT(x)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:,:)
    CALL fftw_execute_dft(pFFT,x,x)
  END SUBROUTINE FFT

  SUBROUTINE IFT(x)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:,:)
    CALL fftw_execute_dft(pIFT,x,x)
  END SUBROUTINE IFT

  SUBROUTINE RFT(x,rx) 
    complex(C_DOUBLE_COMPLEX) :: x(:,:)
    real(C_DOUBLE) :: rx(:,:)
    CALL fftw_execute_dft_c2r(pRFT,x,rx)
  END SUBROUTINE RFT

  SUBROUTINE HFT(rx,x) 
    complex(C_DOUBLE_COMPLEX) :: x(:,:)
    real(C_DOUBLE) :: rx(:,:)
    CALL fftw_execute_dft_r2c(pHFT,rx,x)
  END SUBROUTINE HFT

  SUBROUTINE free_fft()
    CALL fftw_free(ptrVVdx); CALL fftw_free(ptrVVdz);
  END SUBROUTINE free_fft


END MODULE ffts
