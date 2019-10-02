!-----------------------------------------------------------------------
!> @brief
!> Numerical flux is a combination of two fluxes on the interface.
!! Based on the interface is an element interface or a boundary interface,
!! the numerical flux is computed by either two correspond fluxes on the 
!! element interfaces, or element interface and boundary condition (we 
!! call it external state).
!! On each direction, i.e. x and y, each element has two interfaces. 
!! And the target element's neighbor is decided by its dual coordinate in
!! i-j plane.   
!-----------------------------------------------------------------------

MODULE NUMERICAL_FLUX

USE MPI

IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
SUBROUTINE NUMERICAL_FLUX_X(N1, M1)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N1   ! POLY ORDER IN X
    INTEGER, INTENT(IN) :: M1   ! POLY ORDER IN Y
    
    INTEGER :: J
    
    DO J=0, M1
        
    
    ENDDO

END SUBROUTIEN NUMERICAL_FLUX_X

END MODULE NUMERICAL_FLUX
