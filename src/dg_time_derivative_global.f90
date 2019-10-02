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
    
    ! X DIRECTION-------------------------------------------------------
    
    !FIRST GET FLUX ON THE INFERFACES GLOBALLY
    ALLOCATE(SOLUTION_INT_L(0:MMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1))
    ALLOCATE(SOLUTION_INT_R(0:MMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1))
    
    SOLUTION_INT_L = 0.0D0; SOLUTION_INT_R = 0.0D0
    
    ! THEN INTERPOLATES TO BOUNDARY
    DO K = 0, NUM_OF_ELEMENT-1
    
        CALL POLY_LEVEL_TO_ORDER(N, PLEVEL_X(K), PORDER_X)
        CALL POLY_LEVEL_TO_ORDER(M, PLEVEL_Y(K), PORDER_Y)
        
        CALL CONSTRUCT_INTERFACES_X(PORDER_X, PORDER_Y, NUM_OF_EQUATION, &
                                    SOLUTION(0:PORDER_X, 0:PORDER_Y, NUM_OF_EQUATION, K), &
                                    LAGRANGE_LEFT_T(0:PORDER_X, PLEVEL_X(K)), &
                                    LAGRANGE_RIGHT_T(0:PORDER_X, PLEVEL_X(K)), &
                                    SOLUTION_INT_L(0:PORDER_Y, NUM_OF_EQUATION, K), &
                                    SOLUTION_INT_R(0:PORDER_Y, NUM_OF_EQUATION, K) )
        
    
    ENDDO
    
    ! NEXT STEP COMPUTE NUMERICAL FLUXES
    
    
    !-------------------------------------------------------------------
    
    ! Y ----------------------------------------------------------------
    DO K = 0, NUM_OF_ELEMENT-1
        
    
    ENDDO
    !-------------------------------------------------------------------


END SUBROUTINE DG_TIME_DER_COMBINE

END MODULE DG_TIME_DER_ALL
