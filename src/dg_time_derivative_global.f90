!-----------------------------------------------------------------------
!> @brief
!> Compute the dg time direvative globally.
!-----------------------------------------------------------------------

MODULE DG_TIME_DER_ALL

USE MPI
USE NODAL_2D_STORAGE
USE PARAM, ONLY: NUM_OF_EQUATION, NMAX, MMAX, N, M
USE POLY_LEVEL_AND_ORDER
USE INTERFACES_CONSTRUCT
USE NUMERICAL_FLUX
USE A_TIMES_SPATIAL_DER

IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------
!> Calculate the time direvative for current time step,
!! element by element globally. 
!-----------------------------------------------------------------------
SUBROUTINE DG_TIME_DER_COMBINE(T)

    IMPLICIT NONE
    
    INTEGER :: K
    INTEGER :: PORDER_X, PORDER_Y   ! POLY ORDER
    
    DOUBLE PRECISION, INTENT(IN) :: T   !< CURRENT TIME STEP
    
    ! X DIRECTION=======================================================
    
    !FIRST GET FLUX ON THE INFERFACES GLOBALLY--------------------------
    ALLOCATE(SOLUTION_INT_L(0:MMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1))
    ALLOCATE(SOLUTION_INT_R(0:MMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1))
    
    SOLUTION_INT_L = 0.0D0; SOLUTION_INT_R = 0.0D0
    
    
    DO K = 0, NUM_OF_ELEMENT-1
    
        CALL POLY_LEVEL_TO_ORDER(N, PLEVEL_X(K), PORDER_X)
        CALL POLY_LEVEL_TO_ORDER(M, PLEVEL_Y(K), PORDER_Y)
        
        CALL CONSTRUCT_INTERFACES_X(PORDER_X, PORDER_Y, NUM_OF_EQUATION, &
                                    SOLUTION(0:PORDER_X, 0:PORDER_Y, :, K), &
                                    LAGRANGE_LEFT_T(0:PORDER_X, PLEVEL_X(K)), &
                                    LAGRANGE_RIGHT_T(0:PORDER_X, PLEVEL_X(K)), &
                                    SOLUTION_INT_L(0:PORDER_Y, :, K), &
                                    SOLUTION_INT_R(0:PORDER_Y, :, K) )
        
    
    ENDDO
    !-------------------------------------------------------------------
    
    ! NEXT STEP COMPUTE NUMERICAL FLUXES--------------------------------
    ALLOCATE(NFLUX_X_L(0:MMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1))
    ALLOCATE(NFLUX_X_R(0:MMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1))
    
    NFLUX_X_L = 0.0D0; NFLUX_X_R = 0.0D0
    
    
    DO K = 0, NUM_OF_ELEMENT-1
        CALL NUMERICAL_FLUX_X(K, T)
    ENDDO
    !-------------------------------------------------------------------
    
    ! SPACIAL DERIVATIVE------------------------------------------------
    ALLOCATE(FLUX_X(0:NMAX, 0:MMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1))
    
    FLUX_X = 0.0D0
    
    ALLOCATE(FLUX_DER_X(0:NMAX, 0:MMAX, NUM_OF_EQUATION, 0: NUM_OF_ELEMENT-1))
    
    FLUX_DER_X = 0.0D0
    
    ALLOCATE(SOLUTION_TIME_DER(0:NMAX, 0:MMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1))
    
    SOLUTION_TIME_DER = 0.0D0
    
    DO K = 0, NUM_OF_ELEMENT-1
        CALL A_TIMES_SPATIAL_DERIRATIVE_X(K)
    
    ENDDO
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    DEALLOCATE(SOLUTION_INT_L, SOLUTION_INT_R)
    DEALLOCATE(NFLUX_X_L, NFLUX_X_R)
    DEALLOCATE(FLUX_X, FLUX_DER_X)
    !-------------------------------------------------------------------
    
    !===================================================================
    
    
    ! Y ================================================================
    !FIRST GET FLUX ON THE INFERFACES GLOBALLY--------------------------
    ALLOCATE(SOLUTION_INT_L(0:NMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1))
    ALLOCATE(SOLUTION_INT_R(0:NMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1))
    
    SOLUTION_INT_L = 0.0D0; SOLUTION_INT_R = 0.0D0
    
    
    DO K = 0, NUM_OF_ELEMENT-1
    
        CALL POLY_LEVEL_TO_ORDER(N, PLEVEL_X(K), PORDER_X)
        CALL POLY_LEVEL_TO_ORDER(M, PLEVEL_Y(K), PORDER_Y)
        
        CALL CONSTRUCT_INTERFACES_Y(PORDER_X, PORDER_Y, NUM_OF_EQUATION, &
                                    SOLUTION(0:PORDER_X, 0:PORDER_Y, :, K), &
                                    LAGRANGE_DOWN_T(0:PORDER_Y, PLEVEL_Y(K)), &
                                    LAGRANGE_UP_T(0:PORDER_Y, PLEVEL_Y(K)), &
                                    SOLUTION_INT_L(0:PORDER_X, :, K), &
                                    SOLUTION_INT_R(0:PORDER_X, :, K) )
        
    
    ENDDO
    !-------------------------------------------------------------------
    
    ! NEXT STEP COMPUTE NUMERICAL FLUXES--------------------------------
    ALLOCATE(NFLUX_Y_D(0:NMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1))
    ALLOCATE(NFLUX_Y_U(0:NMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1))
    
    NFLUX_Y_D = 0.0D0; NFLUX_Y_U = 0.0D0
    
    
    DO K = 0, NUM_OF_ELEMENT-1
        CALL NUMERICAL_FLUX_Y(K, T)
    ENDDO
    !-------------------------------------------------------------------
    
    ! SPACIAL DERIVATIVE------------------------------------------------
    ALLOCATE(FLUX_Y(0:NMAX, 0:MMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1))
    
    FLUX_Y = 0.0D0
    
    ALLOCATE(FLUX_DER_Y(0:NMAX, 0:MMAX, NUM_OF_EQUATION, 0: NUM_OF_ELEMENT-1))
    
    FLUX_DER_Y = 0.0D0
    
    DO K = 0, NUM_OF_ELEMENT-1
        CALL A_TIMES_SPATIAL_DERIRATIVE_Y(K)
    
    ENDDO
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    DEALLOCATE(SOLUTION_INT_L, SOLUTION_INT_R)
    DEALLOCATE(NFLUX_Y_D, NFLUX_Y_U)
    DEALLOCATE(FLUX_Y, FLUX_DER_Y)
    !-------------------------------------------------------------------
    
    !===================================================================


END SUBROUTINE DG_TIME_DER_COMBINE

END MODULE DG_TIME_DER_ALL
