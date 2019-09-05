!-----------------------------------------------------------------------
!> @brief
!> Riemann solver. Solve the numerical flux at the boundary
!-----------------------------------------------------------------------

MODULE RIEMANN_SOLVER

USE MPI
USE PARAM, ONLY: NUM_OF_EQUATION, C

IMPLICIT NONE

CONTAINS

SUBROUTINE RIEMANN(Q_L, Q_R, DIRECTION, N_FLUX)
!-----------------------------------------------------------------------
! ALGORITHM 88
! THE DERIVATION OF THE NUMERICAL FLUX IS KNOWN AS RIEMANN PROBLEM
!-----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: DIRECTION(2) !< WAVE DIRECTION, (NX, NY)

    DOUBLE PRECISION :: Q_L(NUM_OF_EQUATION)    !< SOLUTION AT LEFT
    DOUBLE PRECISION :: Q_R(NUM_OF_EQUATION)    !< SOLUTION AT RIGHT
    
    DOUBLE PRECISION :: P_L, P_R    ! PRESSURE
    DOUBLE PRECISION :: U_L, U_R    ! VELOCITY IN X DIRECTION
    DOUBLE PRECISION :: V_L, V_R    ! VELOCITY IN Y DIRECTION
    DOUBLE PRECISION :: W_L, W_R    ! WAVE COMPONENT
    
    DOUBLE PRECISION :: N_FLUX(NUM_OF_EQUATION) !< NUMERICAL FLUX
    
    P_L = Q_L(1); U_L = Q_L(2); V_L = Q_L(3);
    P_R = Q_R(1); U_R = Q_R(2); V_R = Q_R(3);
    
    W_L =  P_L + C * (DIRECTION(1) * U_L + DIRECTION(2) * V_L)
    W_R =  P_R - C * (DIRECTION(1) * U_R + DIRECTION(2) * V_R)
    
    N_FLUX(1) = C * (W_L - W_R) / 2
    N_FLUX(2) = DIRECTION(1) * (W_L + W_R) / 2
    N_FLUX(3) = DIRECTION(2) * (W_L + W_R) / 2
    

END SUBROUTINE RIEMANN

END MODULE RIEMANN_SOLVER
