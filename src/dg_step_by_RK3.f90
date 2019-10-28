!-----------------------------------------------------------------------
!> @brief
!> The integration in time by using low storage third order Runge-Kutta. 
!! The 3rd order Runge-Kutta methond is an explicit time time integration 
!!The 3rd order Runge-Kutta methond is an explicit time time integration 
!! method. So there is a time step limitation, which depends on the
!! eigenvalues of the derivative matrix. 
!-----------------------------------------------------------------------


MODULE TIME_STEP_BY_RK

USE MPI
USE PARAM, ONLY: NMAX, MMAX, NUM_OF_EQUATION, N, M
USE NODAL_2D_STORAGE
USE DG_TIME_DER_ALL
USE POLY_LEVEL_AND_ORDER
USE LOCAL_STORAGE

IMPLICIT NONE

CONTAINS

! ----------------------------------------------------------------------
! THRID ORDER RUNGE-KUTTA INTEGRATION 
! INPUT: TN: CURRENT TIME; TIME_DER: DG TIME DERIVATIVE
! ALGORITHM 62
!-----------------------------------------------------------------------
SUBROUTINE DG_STEP_BY_RK3(TN, DELTA_T)

    IMPLICIT NONE
    
    INTEGER :: I, J, K, L, ELEM_K
    
    INTEGER :: N1, M1
    
    DOUBLE PRECISION :: TN  !< CURRENT TIME
    DOUBLE PRECISION :: T
    
    DOUBLE PRECISION :: G(0:NMAX, 0:MMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1)  ! EQUIVALENT TO K
    DOUBLE PRECISION :: AM(3), BM(3), GM(3) ! COEFFICIENT

    DOUBLE PRECISION :: DELTA_T !< TIME STEP SIZE

    
    AM = (/0.0D0, (-5.0D0/9.0D0), (-153D0/128.0D0)/)
    BM = (/0.0D0, (1.0D0/3.0D0), (3.0D0/4.0D0)/)
    GM = (/(1.0D0/3.0D0), (15.0D0/16.0D0), (8.0D0/15.0D0)/)
    
    G = 0.0D0
    ! ------------------------------------------------------------------
    DO K=1,3
    
        T=TN+BM(K)*DELTA_T

        ! GET TIME DERIVATIVE AT CURRENT TIME POINT
        CALL DG_TIME_DER_COMBINE(T)
        
        
        DO ELEM_K = 0, NUM_OF_ELEMENT-1
        
            CALL POLY_LEVEL_TO_ORDER(N, PLEVEL_X(ELEM_K), N1)
            CALL POLY_LEVEL_TO_ORDER(M, PLEVEL_Y(ELEM_K), M1)
            
        
            DO L=1, NUM_OF_EQUATION
                DO J=0,M1
                    DO I=0, N1
                        G(I, J, L, ELEM_K) = AM(K)*G(I, J, L, ELEM_K) & 
                                            + SOLUTION_TIME_DER(I, J, L, ELEM_K)
                                            
                        SOLUTION(I, J, L, ELEM_K) = SOLUTION(I, J, L, ELEM_K) &
                                            + GM(K)*DELTA_T*G(I, J, L, ELEM_K)
                    ENDDO
                ENDDO
            ENDDO
        ENDDO

        ! DEALLOCATE
        DEALLOCATE(SOLUTION_TIME_DER)
    ENDDO
    !-------------------------------------------------------------------


END SUBROUTINE DG_STEP_BY_RK3

END MODULE TIME_STEP_BY_RK
