!============================================!
!                                            !
!              Type definitions              !
!                  for the                   !
!     Direct Numerical Simulation (DNS)      !
!        of a turbulent channel flow         !
!                                            !
!============================================!
!
! Author: M.Sc. Davide Gatti
! Date  : 28/Jul/2015
!

MODULE typedef

  USE, intrinsic :: iso_c_binding
  IMPLICIT NONE

  TYPE :: RHSTYPE
    complex(C_DOUBLE_COMPLEX) :: eta,d2v
  END TYPE RHSTYPE

  TYPE :: Di
    real(C_DOUBLE), dimension(-2:2) :: d0,d1,d2,d4
  END TYPE Di

  TYPE :: VELOCITY
    complex(C_DOUBLE_COMPLEX) :: u,v,w
  END TYPE VELOCITY

  TYPE :: MOMFLUX
    complex(C_DOUBLE_COMPLEX) :: uu,vv,ww,uv,vw,uw
  END TYPE MOMFLUX

  TYPE :: REALVELOCITY
    real(C_DOUBLE) :: u,v,w
  END TYPE REALVELOCITY

  TYPE :: REALMOMFLUX
    real(C_DOUBLE) :: uu,vv,ww,uv,vw,uw
  END TYPE REALMOMFLUX

END MODULE typedef
