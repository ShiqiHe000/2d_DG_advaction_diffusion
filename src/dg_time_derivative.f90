!-----------------------------------------------------------------------
!> @brief
!> Time derivative in 2d for DG approximation
!-----------------------------------------------------------------------

MODULE DG_TIME_DERIVATIVE

USE MPI
USE PARAM
USE RIEMANN_SOLVER
USE NODAL_2D_STORAGE
USE RIEMANN_SOLVER
USE FLUX_VECTORS
USE EXTERNAL_STATE
USE SPACTIAL_DERIVATIVE
USE CENTRAL_FLUX
USE LAX_FRIEDRICHES


IMPLICIT NONE

CONTAINS

SUBROUTINE DG_TIME_DER(T)
!-----------------------------------------------------------------------
! Algorithm 93
!-----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: I, J, S
    
    DOUBLE PRECISION :: X, Y    ! AIMED POINT
    DOUBLE PRECISION, INTENT(IN) :: T !< CURRENT TIME
    
    ALLOCATE(SOLUTION_INT_L(NUM_OF_EQUATION), SOLUTION_INT_R(NUM_OF_EQUATION))
    ALLOCATE(SOLUTION_EXT_L(NUM_OF_EQUATION), SOLUTION_EXT_R(NUM_OF_EQUATION))
    ALLOCATE(NFLUX_X_L(NUM_OF_EQUATION), NFLUX_X_R(NUM_OF_EQUATION))
    ALLOCATE(FLUX_X(0:N, NUM_OF_EQUATION))
    ALLOCATE(FLUX_DER_X(0:N, NUM_OF_EQUATION))
    ALLOCATE(SOLUTION_TIME_DER(0:N, 0:M, NUM_OF_EQUATION))
    
    SOLUTION_EXT_L = 0.0D0; SOLUTION_EXT_R = 0.0D0
    SOLUTION_INT_L = 0.0D0; SOLUTION_INT_R = 0.0D0
    NFLUX_X_L = 0.0D0; NFLUX_X_R = 0.0D0
    FLUX_X = 0.0D0
    FLUX_DER_X = 0.0D0
    SOLUTION_TIME_DER = 0.0D0

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

        ! IMPOSE EXACT SOLUTION ON THE BOUNDARIES-----------------------
!        CALL EXTERNAL_STATE_GAUSSIAN_EXACT(NUM_OF_EQUATION, &
!                                        SOLUTION_EXT_L, T, -1.0D0, Y, &
!                                        DELTA_X, DELTA_Y)
!        CALL EXTERNAL_STATE_GAUSSIAN_EXACT(NUM_OF_EQUATION, &
!                                        SOLUTION_EXT_R, T,  1.0D0, Y, &
!                                        DELTA_X, DELTA_Y)
        CALL EXTERNAL_SINU(NUM_OF_EQUATION, SOLUTION_EXT_L, T, &
                            X_L, Y, DELTA_X, DELTA_Y)
        CALL EXTERNAL_SINU(NUM_OF_EQUATION, SOLUTION_EXT_R, T, &
                            X_R, Y, DELTA_X, DELTA_Y)
        !---------------------------------------------------------------
    
!        CALL RIEMANN(SOLUTION_EXT_L, SOLUTION_INT_L, NFLUX_X_L, (/1.0D0, 0.0D0/), -1.0D0)
!        CALL RIEMANN(SOLUTION_INT_R, SOLUTION_EXT_R, NFLUX_X_R, (/1.0D0, 0.0D0/),  1.0D0)

        CALL RIEMANN_X(SOLUTION_EXT_L, SOLUTION_INT_L, NFLUX_X_L, -1.0D0)
        CALL RIEMANN_X(SOLUTION_INT_R, SOLUTION_EXT_R, NFLUX_X_R,  1.0D0)

        ! CALCULATE THE FLUX
        DO I=0, N
            CALL XFLUX(SOLUTION(I, J, :), FLUX_X(I, :))

        ENDDO

        ! COMPUTE THE FLUX DERIVATIVE
        CALL DG_SPATIAL_DERIVATIVE(N, NFLUX_X_L, NFLUX_X_R, &
                                    FLUX_X, FLUX_DER_X, M_FIRST_DER_X, &
                                    LAGRANGE_LEFT, LAGRANGE_RIGHT, &
                                    GL_W_X)
        DO S=1, NUM_OF_EQUATION
            DO I=0, N
                SOLUTION_TIME_DER(I, J, S) = - 2.0D0 * FLUX_DER_X(I, S) / DELTA_X
            ENDDO
        
        ENDDO

    ENDDO

    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    DEALLOCATE(NFLUX_X_L, NFLUX_X_R)
    DEALLOCATE(FLUX_X)
    DEALLOCATE(FLUX_DER_X)
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    ALLOCATE(NFLUX_Y_D(NUM_OF_EQUATION), NFLUX_Y_U(NUM_OF_EQUATION))
    ALLOCATE(FLUX_Y(0:M, NUM_OF_EQUATION))
    ALLOCATE(FLUX_DER_Y(0:M, NUM_OF_EQUATION))
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    NFLUX_Y_D = 0.0D0; NFLUX_Y_U = 0.0D0
    FLUX_Y = 0.0D0 
    FLUX_DER_Y = 0.0D0
    SOLUTION_INT_L = 0.0D0; SOLUTION_INT_R = 0.0D0
    SOLUTION_EXT_L = 0.0D0; SOLUTION_EXT_R = 0.0D0
    !-------------------------------------------------------------------
          
    ! Y DIRECTION-------------------------------------------------------
    DO I=0, N
        
        X = GL_POINT_X(I)
        
        DO S=1, NUM_OF_EQUATION
            CALL INTERPOLATE_TO_BOUNDARY(M, SOLUTION(I, :, S), &
                                        LAGRANGE_DOWN, SOLUTION_INT_L(S))
                                        
                                        
            CALL INTERPOLATE_TO_BOUNDARY(M, SOLUTION(I, :, S), &
                                        LAGRANGE_UP, SOLUTION_INT_R(S))
        
        ENDDO

!        CALL EXTERNAL_STATE_GAUSSIAN_EXACT(NUM_OF_EQUATION, &
!                                        SOLUTION_EXT_L, T, X, -1.0D0, &
!                                        DELTA_X, DELTA_Y)
!        CALL EXTERNAL_STATE_GAUSSIAN_EXACT(NUM_OF_EQUATION, &
!                                        SOLUTION_EXT_R, T, X, 1.0D0, &
!                                        DELTA_X, DELTA_Y)
        CALL EXTERNAL_SINU(NUM_OF_EQUATION, SOLUTION_EXT_L, T, &
                            X, Y_L, DELTA_X, DELTA_Y)
        CALL EXTERNAL_SINU(NUM_OF_EQUATION, SOLUTION_EXT_R, T, &
                            X, Y_R, DELTA_X, DELTA_Y)
!        CALL RIEMANN(SOLUTION_EXT_L, SOLUTION_INT_L, NFLUX_Y_D, (/0.0D0, 1.0D0/), -1.0D0)
!        CALL RIEMANN(SOLUTION_INT_R, SOLUTION_EXT_R, NFLUX_Y_U, (/0.0D0, 1.0D0/),  1.0D0)

        CALL RIEMANN_Y(SOLUTION_EXT_L, SOLUTION_INT_L, NFLUX_Y_D, -1.0D0)
        CALL RIEMANN_Y(SOLUTION_INT_R, SOLUTION_EXT_R, NFLUX_Y_U,  1.0D0)

        DO J=0, M
            CALL  YFLUX(SOLUTION(I, J, :), FLUX_Y(J, :))
        
        ENDDO
   
        CALL DG_SPATIAL_DERIVATIVE(M, NFLUX_Y_D, NFLUX_Y_U, &
                                    FLUX_Y, FLUX_DER_Y, M_FIRST_DER_Y, &
                                    LAGRANGE_DOWN, LAGRANGE_UP, &
                                    GL_W_Y)
                            

        DO J=0, M
            DO S=1, NUM_OF_EQUATION
                SOLUTION_TIME_DER(I, J, S) = SOLUTION_TIME_DER(I, J, S) &
                                            - 2.0D0 *FLUX_DER_Y(J, S) / DELTA_Y
            ENDDO
        ENDDO
        
    ENDDO
    !-------------------------------------------------------------------

   !-------------------------------------------------------------------
    DEALLOCATE(NFLUX_Y_D, NFLUX_Y_U)
    DEALLOCATE(FLUX_Y)
    DEALLOCATE(FLUX_DER_Y)
    DEALLOCATE(SOLUTION_INT_L, SOLUTION_INT_R)
    DEALLOCATE(SOLUTION_EXT_L, SOLUTION_EXT_R)
    !-------------------------------------------------------------------

END SUBROUTINE DG_TIME_DER


END MODULE DG_TIME_DERIVATIVE
