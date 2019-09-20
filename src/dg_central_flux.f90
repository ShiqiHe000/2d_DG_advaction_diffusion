!-----------------------------------------------------------------------
!> @brief
!> The central flux is the average of the solutions on the two side of 
!! the interface 
!-----------------------------------------------------------------------

MODULE CENTRAL_FLUX

USE MPI
USE PARAM, ONLY: NUM_OF_EQUATION, C

IMPLICIT NONE


CONTAINS 

SUBROUTINE CENTRAL_NUMERICAL_FLUX(Q_L, Q_R, N_FLUX, NORMAL)

    IMPLICIT NONE
    
    INTEGER :: I
    
    DOUBLE PRECISION :: Q_L(NUM_OF_EQUATION)    !< SOLUTION AT LEFT
    DOUBLE PRECISION :: Q_R(NUM_OF_EQUATION)    !< SOLUTION AT RIGHT
    
    DOUBLE PRECISION :: N_FLUX(NUM_OF_EQUATION) !< NUMERICAL FLUX
    
    DOUBLE PRECISION :: NORMAL  !< NORMAL VECTOR, POINTING OUTWARDS OF THE BOUNDARY
    
    DO I=1, NUM_OF_EQUATION
        N_FLUX(I) = 0.5D0 * (Q_L(I) + Q_R(I)) * NORMAL
    
    ENDDO

END SUBROUTINE CENTRAL_NUMERICAL_FLUX



END MODULE CENTRAL_FLUX
