!-----------------------------------------------------------------------
!> brief
!> Driver for DG approxiation. 
!! First, get DG basis parameters, such as collocation points and weights.
!! Then marches by each time step. Using explicit 3rd order Runge-Kutta 
!! method. 
!-----------------------------------------------------------------------

MODULE ADVECTION_DIFFUSION_DRIVER

USE MPI
USE PARAM, ONLY: N, M, T_TOTAL, NT, NUM_OF_EQUATION, NMAX, MMAX
USE DG_2D_CONSTRUCTOR
USE TIME_STEP_BY_RK
USE USER_DEFINE
USE NODAL_2D_STORAGE, ONLY: NUM_OF_ELEMENT
USE POLY_LEVEL_AND_ORDER

IMPLICIT NONE

CONTAINS

SUBROUTINE DRIVER_FOR_DG_APPROXIMATION
!-----------------------------------------------------------------------
! ALGORITHM 51
!-----------------------------------------------------------------------

    IMPLICIT NONE
    
    INTEGER :: K
    INTEGER :: N_NOW, M_NOW !< CURRENT POLY ORDER
    
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
    ALLOCATE(SOLUTION(0:NMAX, 0:MMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1))
    SOLUTION = 0.0D0
    !-------------------------------------------------------------------
    
    ! INITIALZIE SOLUTION-----------------------------------------------
    DO K = 0, NUM_OF_ELEMENT-1
        
        CALL POLY_LEVEL_TO_ORDER(N, PLEVEL_X(K), N_NOW) ! X DIRECTION POLY ORDER
        CALL POLY_LEVEL_TO_ORDER(M, PLEVEL_Y(K), M_NOW) ! Y DIRECTION POLY ORDER
        
        CALL INITIAL_CONDITION_GAUSSIAN(N_NOW, M_NOW, NUM_OF_EQUATION, &
                                        SOLUTION(0:N_NOW, 0:M_NOW, NUM_OF_EQUATION, K), &
                                        GL_POINT_X_T(0:N_NOW, PLEVEL_X(K)), &
                                        GL_POINT_Y_T(0:M_NOW, PLEVEL_Y(K)), &
                                        X_HILBERT(1, K), Y_HILBERT(1, K), &
                                        DELTA_X(K), DELTA_Y(K))
                                      
    ENDDO
!    CALL INITIAL_SINUSOIDAL(N, M, NUM_OF_EQUATION, SOLUTION, GL_POINT_X, GL_POINT_Y)
    !-------------------------------------------------------------------

    ! TIME MARCHES ON---------------------------------------------------
    DO K = 0, NT-1
        CALL DG_STEP_BY_RK3(TN, DELTA_T)
        TN = (K+1) * DELTA_T
    ENDDO
    !-------------------------------------------------------------------
    

END SUBROUTINE DRIVER_FOR_DG_APPROXIMATION

END MODULE ADVECTION_DIFFUSION_DRIVER
