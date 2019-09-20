!-----------------------------------------------------------------------
!> @brief
!> Riemann solver. Solve the numerical flux at the boundary
!-----------------------------------------------------------------------

MODULE RIEMANN_SOLVER

USE MPI
USE PARAM, ONLY: NUM_OF_EQUATION, C

IMPLICIT NONE

CONTAINS

SUBROUTINE RIEMANN(Q_L, Q_R, N_FLUX, UNIT_VECTOR, OUTWARD)
!-----------------------------------------------------------------------
! ALGORITHM 88
! THE DERIVATION OF THE NUMERICAL FLUX IS KNOWN AS RIEMANN PROBLEM
!-----------------------------------------------------------------------

    IMPLICIT NONE

    DOUBLE PRECISION :: Q_L(NUM_OF_EQUATION)    !< SOLUTION AT LEFT
    DOUBLE PRECISION :: Q_R(NUM_OF_EQUATION)    !< SOLUTION AT RIGHT
    
    DOUBLE PRECISION :: P_L, P_R    ! PRESSURE
    DOUBLE PRECISION :: U_L, U_R    ! VELOCITY IN X DIRECTION
    DOUBLE PRECISION :: V_L, V_R    ! VELOCITY IN Y DIRECTION
    DOUBLE PRECISION :: W_L, W_R    ! WAVE COMPONENT
    
    DOUBLE PRECISION :: N_FLUX(NUM_OF_EQUATION) !< NUMERICAL FLUX
    
    DOUBLE PRECISION :: UNIT_VECTOR(2)  !< COORDINATE AXIS UNIT VECTOR, E.G: X-AXIS (1.0D0, 0)
    DOUBLE PRECISION :: OUTWARD !< BOUNDARY OUTWARD DIRECTION, ON THE LEFT = 1.0D0, ON THE RIGTH = 1.0D0
    
    P_L = Q_L(1); U_L = Q_L(2); V_L = Q_L(3);
    P_R = Q_R(1); U_R = Q_R(2); V_R = Q_R(3);
    
    W_L =  P_L + C * ( UNIT_VECTOR(1) * U_L + UNIT_VECTOR(2) * V_L)
    W_R =  P_R - C * ( UNIT_VECTOR(1) * U_R + UNIT_VECTOR(2) * V_R)
    
    N_FLUX(1) = C * (W_L - W_R) / 2 * OUTWARD
    N_FLUX(2) = UNIT_VECTOR(1) * (W_L + W_R) / 2 * OUTWARD
    N_FLUX(3) = UNIT_VECTOR(2) * (W_L + W_R) / 2 * OUTWARD

END SUBROUTINE RIEMANN

SUBROUTINE RIEMANN_X(Q_L, Q_R, N_FLUX, NORMAL)
!-----------------------------------------------------------------------
! Riemann solver in x direction
! Use upwind fluxes
!-----------------------------------------------------------------------
    IMPLICIT NONE
    
    DOUBLE PRECISION :: Q_L(NUM_OF_EQUATION)    !< SOLUTION AT LEFT
    DOUBLE PRECISION :: Q_R(NUM_OF_EQUATION)    !< SOLUTION AT RIGHT
    
    DOUBLE PRECISION :: N_FLUX(NUM_OF_EQUATION) !< NUMERICAL FLUX
    
    DOUBLE PRECISION :: NORMAL
    
    DOUBLE PRECISION :: P_L, P_R    ! PRESSURE
    DOUBLE PRECISION :: U_L, U_R    ! VELOCITY IN X DIRECTION
    DOUBLE PRECISION :: V_L, V_R    ! VELOCITY IN Y DIRECTION
    DOUBLE PRECISION :: W_L, W_R    ! WAVE COMPONENT
    
    P_L = Q_L(1); U_L = Q_L(2); V_L = Q_L(3);
    P_R = Q_R(1); U_R = Q_R(2); V_R = Q_R(3);
    
    W_L = (P_L + C * U_L) / 2.0D0
    W_R = (P_R - C * U_R) / 2.0D0
    
    N_FLUX(1) = C * (W_L - W_R) * NORMAL
    N_FLUX(2) = (W_L + W_R) * NORMAL
    N_FLUX(3) = 0.0D0
    
END SUBROUTINE RIEMANN_X


SUBROUTINE RIEMANN_Y(Q_L, Q_R, N_FLUX, NORMAL)
!-----------------------------------------------------------------------
! Riemann solver in y direction
! Use upwind fluxes
!-----------------------------------------------------------------------

    IMPLICIT NONE
    
    DOUBLE PRECISION :: Q_L(NUM_OF_EQUATION)    !< SOLUTION AT LEFT
    DOUBLE PRECISION :: Q_R(NUM_OF_EQUATION)    !< SOLUTION AT RIGHT
    
    DOUBLE PRECISION :: N_FLUX(NUM_OF_EQUATION) !< NUMERICAL FLUX
    
    DOUBLE PRECISION :: NORMAL
    
    DOUBLE PRECISION :: P_L, P_R    ! PRESSURE
    DOUBLE PRECISION :: U_L, U_R    ! VELOCITY IN X DIRECTION
    DOUBLE PRECISION :: V_L, V_R    ! VELOCITY IN Y DIRECTION
    DOUBLE PRECISION :: W_L, W_R    ! WAVE COMPONENT

    P_L = Q_L(1); U_L = Q_L(2); V_L = Q_L(3);
    P_R = Q_R(1); U_R = Q_R(2); V_R = Q_R(3);
    
    W_L = (P_L + C * V_L) / 2.0D0
    W_R = (P_R - C * V_R) / 2.0D0
    
    N_FLUX(1) = C * (W_L - W_R) * NORMAL
    N_FLUX(2) = 0.0D0
    N_FLUX(3) = (W_L + W_R) * NORMAL

END SUBROUTINE RIEMANN_Y

END MODULE RIEMANN_SOLVER
