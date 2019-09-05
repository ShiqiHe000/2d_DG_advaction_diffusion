!-----------------------------------------------------------------------
!> @brief
!> The external state reflects a wall/reflection boundary.
!! The reflection boundary is the mirror image of the internal state.
!! Meaning w- wave is created by reflecting the w+ wave at the boundary, 
!! i.e. w- = w+. 
!! The reflection implies that the normal velocity is zero. 
!-----------------------------------------------------------------------

MODULE EXTERNAL_STATE

USE MPI

IMPLICIT NONE

CONTAINS

SUBROUTINE GET_EXTERNAL_STATE(NUM_OF_EQUATION, ALPHA, BETA, Q_INT, DIRECTION, Q_EXT)
!-----------------------------------------------------------------------
! PAGE 212
!-----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NUM_OF_EQUATION    !< NUMBER OF EQUATION
    
    DOUBLE PRECISION, INTENT(IN) :: ALPHA, BETA !< TWO CONSTANT RELATED WITH WAVEVECTOR 
    DOUBLE PRECISION :: Q_INT(NUM_OF_EQUATION)  !< INTERIOR SOLUTIONS ON THE BOUNDARY
    DOUBLE PRECISION :: DIRECTION(2)    !< UNIT VECTOR
    DOUBLE PRECISION :: Q_EXT(NUM_OF_EQUATION)  !< EXTERNAL STATES AT THE BOUNDARY
    
    Q_EXT(1) = Q_INT(1)
    Q_EXT(2) = Q_INT(2) * (BETA**2 - ALPHA**2) - 2.0D0 * ALPHA * BETA * Q_INT(3)
    Q_EXT(3) = -2.0D0 * ALPHA * BETA * Q_INT(2) + (ALPHA**2 - BETA**2) * Q_INT(3) 

END SUBROUTINE GET_EXTERNAL_STATE

END MODULE EXTERNAL_STATE
