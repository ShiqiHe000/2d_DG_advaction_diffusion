!-----------------------------------------------------------------------
!> @brief 
!> This module can transform between polynomial order and level 
!-----------------------------------------------------------------------

MODULE POLY_LEVEL_AND_ORDER

USE MPI

IMPLICIT NONE 

CONTAINS

!-----------------------------------------------------------------------
!> TRANSFORM POLYNOMIAL LEVEL TO POLYNOMIAL ORDER
!-----------------------------------------------------------------------
SUBROUTINE POLY_LEVEL_TO_ORDER(N_MIN, PLEVEL, PORDER)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: N_MIN    !< THE MINIMUM POLYNOMIAL ORDER
    INTEGER, INTENT(IN) :: PLEVEL   !< POLYNOMIAL LEVEL
    INTEGER :: PORDER               !< POLYNOMIAL ORDER
    
    PORDER = N_MIN + 2 * (PLEVEL - 1)


END SUBROUTINE POLY_LEVEL_TO_ORDER

!-----------------------------------------------------------------------
!> TRANSFORM POLYNOMIAL ORDER TO POLYNOMIAL LEVEL
!-----------------------------------------------------------------------
SUBROUTINE POLY_ORDER_TO_LEVEL(N_MIN, PLEVEL, PORDER)
    
    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: N_MIN    !< THE MINIMUM POLYNOMIAL ORDER
    INTEGER :: PLEVEL   !< POLYNOMIAL LEVEL
    INTEGER, INTENT(IN) :: PORDER               !< POLYNOMIAL ORDER
    
    PLEVEL = (PORDER - N_MIN) / 2 + 1

END SUBROUTINE POLY_ORDER_TO_LEVEL

END MODULE POLY_LEVEL_AND_ORDER
