!-----------------------------------------------------------------------
!> brief
!> Driver for DG approxiation. 
!! First, get DG basis parameters, such as collocation points and weights.
!! Then marches by each time step. Using explicit 3rd order Runge-Kutta 
!! method. 
!-----------------------------------------------------------------------

MODULE ADVECTION_DIFFUSION_DRIVER

USE MPI
USE PARAM, ONLY: N, M, T_TOTAL, NT, NUM_OF_EQUATION
USE DG_2D_CONSTRUCTOR
USE TIME_STEP_BY_RK
USE USER_DEFINES

IMPLICIT NONE

CONTAINS

SUBROUTINE DRIVER_FOR_DG_APPROXIMATION
!-----------------------------------------------------------------------
! ALGORITHM 51
!-----------------------------------------------------------------------

    IMPLICIT NONE
    
    INTEGER :: K
    
    DOUBLE PRECISION :: DELTA_T     !< TIME STEP 
    DOUBLE PRECISION :: TN          !< CURRENT TIME
    
    ! CONSTRUCT DG BASIS------------------------------------------------
    CALL CONSTRUCT_BASIS    ! NOW WE HAVE GL POINTS, WEIGHTS, M_FIRST DERIVATIVE MATRICES
    !-------------------------------------------------------------------
    
    ! TIME STEP---------------------------------------------------------
    DELTA_T = T_TOTAL / NT
    !-------------------------------------------------------------------
    
    ! CURRENT TIME------------------------------------------------------
    TN = 0.0D0
    !-------------------------------------------------------------------
    
    ! ALLOCATE----------------------------------------------------------
    ALLOCATE(SOLUTION(0:N, 0:M, NUM_OF_EQUATION))
    !-------------------------------------------------------------------
    
    ! INITIALZIE SOLUTION-----------------------------------------------
    CALL INITIAL_CONDITION(N, M, NUM_OF_EQUATION, SOLUTION, GL_POINT_X, GL_POINT_Y)
    !-------------------------------------------------------------------
    
    ! TIME MARCHES ON---------------------------------------------------
    DO K = 0, NT-1
        CALL DG_STEP_BY_RK3(TN, DELTA_T)
        TN = (K+1) * DELTA_T
    ENDDO
    !-------------------------------------------------------------------
    

END SUBROUTINE DRIVER_FOR_DG_APPROXIMATION

END MODULE ADVECTION_DIFFUSION_DRIVER
