!-----------------------------------------------------------------------
!> brief
!> Driver for DG approxiation. 
!! First, get DG basis parameters, such as collocation points and weights.
!! Then marches by each time step. Using explicit 3rd order Runge-Kutta 
!! method. 
!-----------------------------------------------------------------------

MODULE ADVECTION_DIFFUSION_DRIVER

USE MPI
USE PARAM, ONLY: N, M, T_TOTAL, NT
USE DG_2D_CONSTRUCTOR
USE TIME_STEP_BY_RK

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
    
    ! INITIALZIE SOLUTION-----------------------------------------------
    
    !-------------------------------------------------------------------
    
    DO K = 0, NT-1
        CALL DG_STEP_BY_RK3(TN, DELTA_T)
        TN = (K+1) * DELTA_T
    ENDDO
    

END SUBROUTINE DRIVER_FOR_DG_APPROXIMATION

END MODULE ADVECTION_DIFFUSION_DRIVER
