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
USE PARAM, ONLY: C, K_X, K_Y, Y_L, X_L  ! EXCEPT NUM_OF_EQUATION
USE AFFINE_MAP

IMPLICIT NONE

CONTAINS

SUBROUTINE EXTERNAL_STATE_GAUSSIAN_REFLECT(NUM_OF_EQUATION, ALPHA,&
                                            BETA, Q_INT, Q_EXT)
!-----------------------------------------------------------------------
! PAGE 212
!-----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NUM_OF_EQUATION    !< NUMBER OF EQUATION
    
    DOUBLE PRECISION, INTENT(IN) :: ALPHA, BETA !< TWO CONSTANT RELATED WITH WAVEVECTOR 
    DOUBLE PRECISION :: Q_INT(NUM_OF_EQUATION)  !< INTERIOR SOLUTIONS ON THE BOUNDARY
    DOUBLE PRECISION :: Q_EXT(NUM_OF_EQUATION)  !< EXTERNAL STATES AT THE BOUNDARY
        
    Q_EXT(:) = 0.0D0    

    Q_EXT(1) = Q_INT(1)
    Q_EXT(2) = Q_INT(2) * (BETA**2 - ALPHA**2) - 2.0D0 * ALPHA * BETA * Q_INT(3)
    Q_EXT(3) = -2.0D0 * ALPHA * BETA * Q_INT(2) + (ALPHA**2 - BETA**2) * Q_INT(3) 
    


END SUBROUTINE EXTERNAL_STATE_GAUSSIAN_REFLECT

SUBROUTINE EXTERNAL_STATE_GAUSSIAN_EXACT(NUM_OF_EQUATION, Q_EXT, T, XR, YR, &
                                            DEL_X, DEL_Y)
!-----------------------------------------------------------------------
! IMPLEMENT EXACT SOLUTION ON THE 4 SIDES OF THE BOUNDARIES
!-----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NUM_OF_EQUATION    !< NUMBER OF EQUATION
     
    DOUBLE PRECISION :: Q_EXT(NUM_OF_EQUATION)  !< EXTERNAL STATES AT THE BOUNDARY
    DOUBLE PRECISION, INTENT(IN) :: T   !< CURRENT TIME 
    DOUBLE PRECISION, INTENT(IN) :: XR, YR    !< COORDINATE ON THE BOUNDARY(REFERENCE SPACE)
    DOUBLE PRECISION,INTENT(IN) :: DEL_Y, DEL_X !< ELEMENT SIZE
    
    DOUBLE PRECISION :: INTER   ! INTERMIDIATE VARIABLE
    DOUBLE PRECISION :: X, Y    ! COORDINATES ON THE PHYSICAL SPACE
    
    CALL AFFINE_MAPPING(XR, X, X_L, DEL_X)
    CALL AFFINE_MAPPING(YR, Y, Y_L, DEL_Y)
    
    INTER = DEXP( - (K_X * (X - X0) + &
                            K_Y * (Y - Y0) - C * T)**2 / D**2)
            
    Q_EXT(1) = INTER 
    
    Q_EXT(2) = K_X / C * INTER 
    
    Q_EXT(3) = K_Y / C * INTER 
    
END SUBROUTINE EXTERNAL_STATE_GAUSSIAN_EXACT


SUBROUTINE EXTERNAL_SINU(NUM_OF_EQUATION, Q_EXT, T, XR, YR, DEL_X, DEL_Y)


    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NUM_OF_EQUATION    !< NUMBER OF EQUATION
     
    DOUBLE PRECISION :: Q_EXT(NUM_OF_EQUATION)  !< EXTERNAL STATES AT THE BOUNDARY
    DOUBLE PRECISION, INTENT(IN) :: T   !< CURRENT TIME 
    DOUBLE PRECISION, INTENT(IN) :: XR, YR    !< COORDINATE ON THE BOUNDARY(REFERENCE)
    DOUBLE PRECISION,INTENT(IN) :: DEL_Y, DEL_X !< ELEMENT SIZE
    
    DOUBLE PRECISION :: X, Y    ! COORDINATES ON THE PHYSICAL SPACE
    
    CALL AFFINE_MAPPING(XR, X, X_L, DEL_X)
    CALL AFFINE_MAPPING(YR, Y, Y_L, DEL_Y)
    
    Q_EXT(1) = C * DCOS(X + Y - C*T)
    Q_EXT(2) = DCOS(X + Y - C*T)
    Q_EXT(3) = DCOS(X + Y - C*T)
    

END SUBROUTINE EXTERNAL_SINU

END MODULE EXTERNAL_STATE
