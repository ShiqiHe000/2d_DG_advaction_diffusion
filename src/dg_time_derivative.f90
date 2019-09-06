!-----------------------------------------------------------------------
!> @brief
!> Time derivative in 2d for DG approximation
!-----------------------------------------------------------------------

MODULE DG_TIME_DERIVATIVE

USE MPI
USE PARAM, ONLY: N, M, NUM_OF_EQUATION, K_X, K_Y
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
    ALLOCATE(SOLUTION_EXT_L(NUM_OF_EQUATION), SOLUTION_EXT_R(NUM_OF_EQUATION))
    ALLOCATE(NFLUX_X_L(NUM_OF_EQUATION), NFLUX_X_R(NUM_OF_EQUATION))
    ALLOCATE(FLUX_X(0:N, NUM_OF_EQUATION))
    ALLOCATE(FLUX_DER_X(0:N, NUM_OF_EQUATION))
    ALLOCATE(SOLUTION_TIME_DER(0:N, 0:M, NUM_OF_EQUATION))
    
    ! X DIRECTION-------------------------------------------------------
    DO J=0, M
    
        Y = GL_POINT_Y(J)
        
        ! INTERPOLATE SOLUTION TO THE BOUNDARIES
        DO S=1, NUM_OF_EQUATION
            CALL INTERPOLATE_TO_BOUNDARY(N, SOLUTION(:, J, S), &
                                        LAGRANGE_LEFT, SOLUTION_INT_L(S))
                                        
            CALL INTERPOLATE_TO_BOUNDARY(N, SOLUTION(:, J, S), &
                                        LAGRANGE_RIGHT, SOLUTION_INT_R(S))
        
        ENDDO
        
        
        CALL GET_EXTERNAL_STATE(NUM_OF_EQUATION, K_X, K_Y, &
                                SOLUTION_INT_L, SOLUTION_EXT_L)
                                
        CALL GET_EXTERNAL_STATE(NUM_OF_EQUATION, K_X, K_Y, &
                                SOLUTION_INT_R, SOLUTION_EXT_R)
                                
        CALL RIEMANN(SOLUTION_INT_L, SOLUTION_EXT_L, NFLUX_X_L, K_X, K_Y)
        CALL RIEMANN(SOLUTION_INT_R, SOLUTION_EXT_R, NFLUX_X_R, K_X, K_Y)
        
        ! CALCULATE THE FLUX
        DO I=0, N
            CALL XFLUX(SOLUTION(I, J, :), FLUX_X(I, :))
        
        ENDDO
        
        ! COMPUTE THE FLUX DERIVATIVE
        CALL DG_SPACTIAL_DERIVATIVE(N, NFLUX_X_L, NFLUX_X_R, &
                                    FLUX_X, FLUX_DER_X, M_FIRST_DER_X, &
                                    LAGRANGE_LEFT, LAGRANGE_RIGHT, &
                                    GL_W_X)
                                    
        DO I=0, N
            DO S=1, NUM_OF_EQUATION
                SOLUTION_TIME_DER(I, J, S) = - FLUX_DER_X(I, S)
            ENDDO
        
        ENDDO
    
    ENDDO

END SUBROUTINE DG_TIME_DER


END MODULE DG_TIME_DERIVATIVE
