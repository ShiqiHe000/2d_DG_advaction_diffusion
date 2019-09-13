!-----------------------------------------------------------------------
!> @brief
!> Riemann solver. Solve the numerical flux at the boundary
!-----------------------------------------------------------------------

MODULE RIEMANN_SOLVER

USE MPI
USE PARAM, ONLY: NUM_OF_EQUATION, C

IMPLICIT NONE

CONTAINS

SUBROUTINE RIEMANN(Q_L, Q_R, N_FLUX, ALPHA, BETA)
!-----------------------------------------------------------------------
! ALGORITHM 88
! THE DERIVATION OF THE NUMERICAL FLUX IS KNOWN AS RIEMANN PROBLEM
!-----------------------------------------------------------------------

    IMPLICIT NONE
    
    DOUBLE PRECISION, INTENT(IN) :: ALPHA, BETA !< TWO CONSTANT RELATED WITH WAVEVECTOR 

    DOUBLE PRECISION :: Q_L(NUM_OF_EQUATION)    !< SOLUTION AT LEFT
    DOUBLE PRECISION :: Q_R(NUM_OF_EQUATION)    !< SOLUTION AT RIGHT
    
    DOUBLE PRECISION :: P_L, P_R    ! PRESSURE
    DOUBLE PRECISION :: U_L, U_R    ! VELOCITY IN X DIRECTION
    DOUBLE PRECISION :: V_L, V_R    ! VELOCITY IN Y DIRECTION
    DOUBLE PRECISION :: W_L, W_R    ! WAVE COMPONENT
    
    DOUBLE PRECISION :: N_FLUX(NUM_OF_EQUATION) !< NUMERICAL FLUX
    
    P_L = Q_L(1); U_L = Q_L(2); V_L = Q_L(3);
    P_R = Q_R(1); U_R = Q_R(2); V_R = Q_R(3);
    
    W_L =  P_L + C * (ALPHA * U_L + BETA * V_L)
    W_R =  P_R - C * (ALPHA * U_R + BETA * V_R)
    
    N_FLUX(1) = C * (W_L - W_R) / 2
    N_FLUX(2) = ALPHA * (W_L + W_R) / 2
    N_FLUX(3) = BETA * (W_L + W_R) / 2

!print *, "u_l", u_l, "u_r", u_r
!print *, "v_l", v_l, "v_r", v_r
!!    print *, "w_l", w_l, "w_r", w_r
!    print *, "n_flux", n_flux

END SUBROUTINE RIEMANN

END MODULE RIEMANN_SOLVER
