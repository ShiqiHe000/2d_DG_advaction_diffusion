!-----------------------------------------------------------------------
!< @brief
!< Verify your results!
!-----------------------------------------------------------------------

MODULE VERIFICATION

USE MPI
USE PARAM, ONLY: N, M, NUM_OF_EQUATION, T_TOTAL
USE USER_DEFINES
USE NODAL_2D_STORAGE, ONLY: GL_POINT_X, GL_POINT_Y, SOLUTION

IMPLICIT NONE

CONTAINS

SUBROUTINE GET_ERROR(EXACT, ERROR, L2_NORM)

    IMPLICIT NONE
    
    INTEGER :: I, J, K
    
    DOUBLE PRECISION :: EXACT(0:N, 0:M, NUM_OF_EQUATION)    !< EXACT SOLUTION
    DOUBLE PRECISION :: ERROR(0:N, 0:M, NUM_OF_EQUATION)    !< ERRORS AT COLLOCATION POINTS
    
    DOUBLE PRECISION :: L2_NORM(NUM_OF_EQUATION)    !< WE EVALUATE THE ERRORS AS L2 NORMS
    
    ! GET EXACT SOLUTION------------------------------------------------
    CALL EXACT_SOLUTION(N, M, NUM_OF_EQUATION, GL_POINT_X, GL_POINT_Y, EXACT, T_TOTAL)
    !-------------------------------------------------------------------
    
    ! COMPUET ERROR-----------------------------------------------------
    DO K=1, NUM_OF_EQUATION
        DO J=0, M
            DO I=0, N
                ERROR(I, J, K) = DABS(SOLUTION(I, J, K) - EXACT(I, J, K))
            ENDDO
        
        ENDDO
    ENDDO
    !-------------------------------------------------------------------
    
    ! L2 NORMS----------------------------------------------------------
    DO K=1, NUM_OF_EQUATION
        L2_NORM(K) = NORM2(ERROR(:, :, K))
    ENDDO
    !-------------------------------------------------------------------

END SUBROUTINE GET_ERROR


END MODULE VERIFICATION
