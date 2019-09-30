!-----------------------------------------------------------------------
!> @brief
!> Generate element k's dual graph coorindate.
!! Constraint: physical domain should be in first quadrant in Cartesian
!! plane.
!-----------------------------------------------------------------------

MODULE GET_DUAL_COORD

USE MPI

IMPLICIT NONE 

CONTAINS

SUBROUTINE GET_DUAL_COORD_2D(X, Y, DELTA_X, DELTA_Y, COORD)

    IMPLICIT NONE 
    
    DOUBLE PRECISION, INTENT(IN) :: X, Y !< ELEMENT K NODE 1 COORD
    
    DOUBLE PRECISION, INTENT(IN) :: DELTA_X, DELTA_Y    !< ELEMENT SIZE
    
    INTEGER :: COORD(2)    !< OUTPUT: ELEMENT DUAL COORD

    COORD(1) = NINT(X / DELTA_X)
    COORD(2) = NINT(Y / DELTA_Y)
    
!        PRINT *, X, DELTA_X
!      PRINT *, X / DELTA_X, COORD(1)

END SUBROUTINE GET_DUAL_COORD_2D


END MODULE GET_DUAL_COORD
