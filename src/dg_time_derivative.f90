!-----------------------------------------------------------------------
!> @brief
!> Time derivative in 2d for DG approximation
!-----------------------------------------------------------------------

MODULE DG_TIME_DERIVATIVE

USE MPI
USE PARAM, ONLY: N, M, NUM_OF_EQUATION
USE RIEMANN_SOLVER
USE NODAL_2D_STORAGE
USE RIEMANN_SOLVER
USE FLUX_VECTORS

IMPLICIT NONE

CONTAINS

SUBROUTINE DG_TIME_DER 
!-----------------------------------------------------------------------
! Algorithm 93
!-----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: I, J, S
    
    DOUBLE PRECISION :: X, Y    ! AIMED POINT
    
    ALLOCATE(SOLUTION(0:N, 0:M, NUM_OF_EQUATION))
    ALLOCATE(SOLUTION_INT_L(NUM_OF_EQUATION), SOLUTION_INT_R(NUM_OF_EQUATION))
    
    ! X DIRECTION-------------------------------------------------------
    DO J=0, M
    
        Y = GL_POINT_Y(J)
        
        DO S=1, NUM_OF_EQUATION
            CALL INTERPOLATE_TO_BOUNDARY(N, SOLUTION(:, J, S), &
                                        LAGRANGE_LEFT, SOLUTION_INT_L(S))
                                        
            CALL INTERPOLATE_TO_BOUNDARY(N, SOLUTION(:, J, S), &
                                        LAGRANGE_RIGHT, SOLUTION_INT_R(S))
        
        ENDDO
        
        
    
    ENDDO

END SUBROUTINE DG_TIME_DER


END MODULE DG_TIME_DERIVATIVE
