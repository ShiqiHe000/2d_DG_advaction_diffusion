!-----------------------------------------------------------------------
!> @brief
!> Flux vectors for the two dimensional wave equation. 
!! The procedure computes the horizontal flux and vertical flux at the 
!! internal grid points. 
!-----------------------------------------------------------------------

MODULE FLUX_VECTORS

USE MPI
USE PARAM, ONLY: NUM_OF_EQUATION, C

IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------
! ALGORITHM 94
!-----------------------------------------------------------------------
SUBROUTINE XFLUX(Q, XF)

    IMPLICIT NONE
    
    DOUBLE PRECISION, INTENT(IN) :: Q(NUM_OF_EQUATION)  !< SOLUTION AT THE AIMED POINT
    DOUBLE PRECISION :: XF(NUM_OF_EQUATION)     !< HORIZONTAL FLUXES
    
    XF(1) = C**2 * Q(2)
    XF(2) = Q(1)
    XF(3) = 0.0D0

    

END SUBROUTINE XFLUX

!-----------------------------------------------------------------------
! ALGORITHM 94
!-----------------------------------------------------------------------
SUBROUTINE YFLUX(Q, YF)

    IMPLICIT NONE
    
    DOUBLE PRECISION, INTENT(IN) :: Q(NUM_OF_EQUATION) !< SOLUTION AT THE AIMED POINT
    DOUBLE PRECISION :: YF(NUM_OF_EQUATION) !< VERTICAL FLUXES
    
    YF(1) = C**2 * Q(3)
    YF(2) = 0.0D0
    YF(3) = Q(1)


END SUBROUTINE YFLUX

END MODULE FLUX_VECTORS
