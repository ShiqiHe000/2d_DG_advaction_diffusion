!-----------------------------------------------------------------------
!> @brief
!> The integration in time by using long storage third order Runge-Kutta. 
!!The 3rd order Runge-Kutta methond is an explicit time time integration 
!! method. So there is a time step limitation, which depends on the
!! eigenvalues of the derivative matrix. 
!-----------------------------------------------------------------------


MODULE TIME_STEP_BY_RK

USE MPI
USE PARAM, ONLY: N, M, NUM_OF_EQUATION
USE DG_TIME_DERIVATIVE
USE NODAL_2D_STORAGE

IMPLICIT NONE

CONTAINS

SUBROUTINE DG_STEP_BY_RK3(TN, DELTA_T)
! ----------------------------------------------------------------------
! THRID ORDER RUNGE-KUTTA INTEGRATION 
! INPUT: TN: CURRENT TIME; TIME_DER: DG TIME DERIVATIVE
! ALGORITHM 62
!-----------------------------------------------------------------------

    IMPLICIT NONE
    
    INTEGER :: I, J, K, L
    
    DOUBLE PRECISION :: TN  !< CURRENT TIME
    DOUBLE PRECISION :: T
    
    DOUBLE PRECISION :: G(0:N, 0:M, NUM_OF_EQUATION)  ! EQUIVALENT TO K
    DOUBLE PRECISION :: AM(3), BM(3), GM(3) ! COEFFICIENT

    DOUBLE PRECISION :: DELTA_T !< TIME STEP SIZE

    
    AM = (/0.0D0, (-5.0D0/9.0D0), (-153D0/128.0D0)/)
    BM = (/0.0D0, (1.0D0/3.0D0), (3.0D0/4.0D0)/)
    GM = (/(1.0D0/3.0D0), (15.0D0/16.0D0), (8.0D0/15.0D0)/)
                                    
    
    ! ------------------------------------------------------------------
    DO K=1,3
    
        T=TN+BM(K)*DELTA_T

        ! GET TIME DERIVATIVE AT CURRENT TIME POINT
        CALL DG_TIME_DER(T)
        
        DO L=1, NUM_OF_EQUATION
            DO J=0,M
                DO I=0, N
                        G(I, J, L) = AM(K)*G(I, J, L) + SOLUTION_TIME_DER(I, J, L)
                        SOLUTION(I, J, L) = SOLUTION(I, J, L) + GM(K)*DELTA_T*G(I, J, L)
                ENDDO
            ENDDO
        ENDDO

        ! DEALLOCATE
        DEALLOCATE(SOLUTION_TIME_DER)
    ENDDO
    !-------------------------------------------------------------------


END SUBROUTINE DG_STEP_BY_RK3

END MODULE TIME_STEP_BY_RK
