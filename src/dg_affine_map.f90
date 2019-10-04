!-----------------------------------------------------------------------
!> @brief
!> We map each element to the reference interval [-1, 1]. which will 
!! serve as our refernece element
!-----------------------------------------------------------------------

MODULE AFFINE_MAP

USE MPI

IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------
! OUTPUT: PHYSICAL LOCATION X
! INPUT: OTHERS
!-----------------------------------------------------------------------
SUBROUTINE AFFINE_MAPPING(XI, X, XK_1, DELTA_X)

    IMPLICIT NONE
    
    DOUBLE PRECISION :: XI  !< REFERENCE LOCATION
    DOUBLE PRECISION :: X   !< PHYSICAL LOCATION, OUTPUT
    DOUBLE PRECISION :: XK_1    !< THE LOCATION OF THE LEFT BOUND OF THIS ELEMENT
    DOUBLE PRECISION :: DELTA_X !< ELEMENT SIZE
    
    X = XK_1 + (XI + 1.0D0) * DELTA_X / 2.0D0


END SUBROUTINE AFFINE_MAPPING

END MODULE AFFINE_MAP
