!-----------------------------------------------------------------------
!< @brief
!< Verify your results!
!-----------------------------------------------------------------------

MODULE VERIFICATION

USE MPI
USE PARAM, ONLY: N, M, NUM_OF_EQUATION, T_TOTAL
USE USER_DEFINES
USE NODAL_2D_STORAGE, ONLY: GL_POINT_X, GL_POINT_Y, SOLUTION
USE WRITE_DATA

IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :, :) :: EXACT  !< EXACT SOLUTION
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :, :) :: ERROR  !< ERRORS AT COLLOCATION POINTS
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: L2_NORM  !< WE EVALUATE THE ERRORS AS L2 NORMS

CONTAINS

SUBROUTINE GET_ERROR
!-----------------------------------------------------------------------
! OUTPUT: EXACT, ERROR, L2_NORM
!-----------------------------------------------------------------------

    IMPLICIT NONE
    
    INTEGER :: I, J, K
    
    !-------------------------------------------------------------------
    ALLOCATE(EXACT(0:N, 0:M, NUM_OF_EQUATION))
    ALLOCATE(ERROR(0:N, 0:M, NUM_OF_EQUATION))
    ALLOCATE(L2_NORM(NUM_OF_EQUATION))
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    EXACT = 0.0D0
    ERROR = 0.0D0
    L2_NORM = 0.0D0
    !-------------------------------------------------------------------
    
    ! GET EXACT SOLUTION------------------------------------------------
!    CALL EXACT_SOLUTION(N, M, NUM_OF_EQUATION, GL_POINT_X, GL_POINT_Y, EXACT, T_TOTAL)
    CALL SIN_EXACT(N, M, NUM_OF_EQUATION, GL_POINT_X, GL_POINT_Y, EXACT, T_TOTAL)
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
    
    CALL WRITE_ERROR(N, M, ERROR(:, :, 2))
    CALL WRITE_RESULTS(N, M, EXACT(:, :, 2), SOLUTION(:, :, 2))
    
    print *, maxval(ERROR(:, 0, 2))
    print *, L2_NORM(2) 
    
END SUBROUTINE GET_ERROR


END MODULE VERIFICATION
