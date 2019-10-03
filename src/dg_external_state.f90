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
USE USER_DEFINES
USE PARAM, ONLY: C
USE USER_DEFINES, ONLY: K_X, K_Y

IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------
! PAGE 212
!-----------------------------------------------------------------------
SUBROUTINE EXTERNAL_STATE_GAUSSIAN_REFLECT(NUM_OF_EQUATION, &
                                            Q_INT, Q_EXT, UNIT_VECTOR)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NUM_OF_EQUATION    !< NUMBER OF EQUATION
    
    DOUBLE PRECISION :: Q_INT(NUM_OF_EQUATION)  !< INTERIOR SOLUTIONS ON THE BOUNDARY
    DOUBLE PRECISION :: Q_EXT(NUM_OF_EQUATION)  !< EXTERNAL STATES AT THE BOUNDARY
          
    DOUBLE PRECISION :: UNIT_VECTOR(2)  !< BOUNDARY NORMAL UNIT VECTOR, DO NOT NEED TO CONSIDER POINTING TO OUTWARD

    Q_EXT(1) = Q_INT(1)
    Q_EXT(2) = - Q_INT(2) * UNIT_VECTOR(1)
    Q_EXT(3) = - Q_INT(3) * UNIT_VECTOR(2)
    


END SUBROUTINE EXTERNAL_STATE_GAUSSIAN_REFLECT

!-----------------------------------------------------------------------
! IMPLEMENT EXACT SOLUTION ON THE 4 SIDES OF THE BOUNDARIES
!-----------------------------------------------------------------------
SUBROUTINE EXTERNAL_STATE_GAUSSIAN_EXACT(NUM_OF_EQUATION, Q_EXT, T, X, Y)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NUM_OF_EQUATION    !< NUMBER OF EQUATION
     
    DOUBLE PRECISION :: Q_EXT(NUM_OF_EQUATION)  !< EXTERNAL STATES AT THE BOUNDARY
    DOUBLE PRECISION, INTENT(IN) :: T   !< CURRENT TIME 
    DOUBLE PRECISION, INTENT(IN) :: X, Y    !< COORDINATE ON THE BOUNDARY
    
    DOUBLE PRECISION :: INTER   ! INTERMIDIATE VARIABLE
    
    INTER = DEXP( - (K_X * (X - X0) + &
                            K_Y * (Y - Y0) - C * T)**2 / D**2)
            
    Q_EXT(1) = INTER 
    
    Q_EXT(2) = K_X / C * INTER 
    
    Q_EXT(3) = K_Y / C * INTER 
    
END SUBROUTINE EXTERNAL_STATE_GAUSSIAN_EXACT

!-----------------------------------------------------------------------
! EXTERNAL STATE FOR SINUSOIDAL CASE
!-----------------------------------------------------------------------
SUBROUTINE EXTERNAL_SINU(NUM_OF_EQUATION, Q_EXT ,X, Y, T)


    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NUM_OF_EQUATION    !< NUMBER OF EQUATION
     
    DOUBLE PRECISION :: Q_EXT(NUM_OF_EQUATION)  !< EXTERNAL STATES AT THE BOUNDARY
    DOUBLE PRECISION, INTENT(IN) :: T   !< CURRENT TIME 
    DOUBLE PRECISION, INTENT(IN) :: X, Y    !< COORDINATE ON THE BOUNDARY
    
    Q_EXT(1) = C * DCOS(X - C*T)
    Q_EXT(2) = DCOS(X - C*T)
    Q_EXT(3) = 0.0D0
    

END SUBROUTINE EXTERNAL_SINU



END MODULE EXTERNAL_STATE
