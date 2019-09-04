!-----------------------------------------------------------------------
!> @brief
!> Compute the first spacial derivative given by the interior points
!! and two boundary solutions. 
!-----------------------------------------------------------------------

MODULE SPACTIAL_DERIVATIVE

USE MPI
USE PARAM, ONLY: NUM_OF_EQUATION
USE BASIS
USE DG_2D_CONSTRUCTOR
USE NODAL_2D_STORAGE

IMPLICIT NONE

!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: FLUX_X_DER    !< FLUX DERIVATIVE
!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: FLUX_Y_DER    !< FLUX DERIVATIVE

CONTAINS

SUBROUTINE DG_SPACTIAL_DERIVATIVE(N_TH, FLUX_LEFT, FLUX_RIGHT, &
                                    FLUX, FLUX_DER, DER, LAG1, LAG2, WEIGHT)
!-----------------------------------------------------------------------
! ALGORITHM 92
! RETURN: FLUX_DER
!-----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: J, S
    INTEGER, INTENT(IN) :: N_TH !< N-TH ORDER POLYNOMIAL
    
    DOUBLE PRECISION :: FLUX_LEFT(NUM_OF_EQUATION)  !< FLUX AT THE END OF BOUNDARY
    DOUBLE PRECISION :: FLUX_RIGHT(NUM_OF_EQUATION) !< FLUX AT THE END OF BOUNDARY
    DOUBLE PRECISION :: FLUX(0:N_TH, NUM_OF_EQUATION)  !< NUMERICAL FLUX
    DOUBLE PRECISION :: FLUX_DER(0:N_TH, NUM_OF_EQUATION)    !< FLUX DERIVATIVE
    DOUBLE PRECISION :: DER(0:N_TH, 0:N_TH) !< FIRST DERIVATIVE MATRIX
    DOUBLE PRECISION :: LAG1(0:N_TH)    !< LAGRANGE INTERPOLATE AT LEFT/BOTTOM
    DOUBLE PRECISION :: LAG2(0:N_TH)    !< LAGRANGE INTERPOLATE AT RIGHT/TOP
    DOUBLE PRECISION :: WEIGHT(0:N_TH)  !< GL WEIGHTS
    
    !-------------------------------------------------------------------
!    ALLOCATE(FLUX_X(0:N, NUM_OF_EQUATION), FLUX_Y(0:M, NUM_OF_EQUATION))
!    ALLOCATE(FLUX_X_DER(0:N, NUM_OF_EQUATION), FLUX_Y_DER(0:M, NUM_OF_EQUATION))
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    DO S=1, NUM_OF_EQUATION
        CALL MATRIX_VECTOR_DERIVATIVE(N_TH, DER, FLUX(:, S), FLUX_DER(:, S))
    ENDDO
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    DO J=0, N_TH
        DO S=1, NUM_OF_EQUATION
            FLUX_DER(J, S) = FLUX_DER(J, S) + &
                                (FLUX_RIGHT(S) * LAG1(J) + &
                                FLUX_LEFT(S) * LAG2(J)) / WEIGHT(J)
        ENDDO
    ENDDO
    !-------------------------------------------------------------------

END SUBROUTINE DG_SPACTIAL_DERIVATIVE

END MODULE SPACTIAL_DERIVATIVE
