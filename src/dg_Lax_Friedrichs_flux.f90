!-----------------------------------------------------------------------
!> @ brief
!> Lax-Friedrichs flux
!-----------------------------------------------------------------------

MODULE LAX_FRIEDRICHES

USE MPI
USE PARAM, ONLY: NUM_OF_EQUATION, C


IMPLICIT NONE

CONTAINS

SUBROUTINE LAX_FRIEDRICHES_FLUX(Q_L, Q_R, N_FLUX, ALPHA, NORMAL)

    IMPLICIT NONE
    
    INTEGER :: I
    
    DOUBLE PRECISION :: Q_L(NUM_OF_EQUATION)    !< SOLUTION AT LEFT
    DOUBLE PRECISION :: Q_R(NUM_OF_EQUATION)    !< SOLUTION AT RIGHT
    
    DOUBLE PRECISION :: N_FLUX(NUM_OF_EQUATION) !< NUMERICAL FLUX
    
    DOUBLE PRECISION :: ALPHA   !< PARAMETER: 0 <= ALPHA <= 1.0
    DOUBLE PRECISION :: NORMAL  !< NORMAL VECTOR, POINTING OUTWARDS OF THE BOUNDARY
    
    DO I=1, NUM_OF_EQUATION
        N_FLUX(I) = 0.5D0 * C * NORMAL * (Q_L(I) + Q_R(I) + &
                        NORMAL * (1.0D0 - ALPHA) * (Q_L(I) - Q_R(I)))
    
    ENDDO

    
    

END SUBROUTINE LAX_FRIEDRICHES_FLUX


END MODULE LAX_FRIEDRICHES
