!-----------------------------------------------------------------------
!> @brief
!> Verify your results!
!-----------------------------------------------------------------------

MODULE VERIFICATION

USE MPI
USE PARAM, ONLY: N, M, NUM_OF_EQUATION, T_TOTAL, NMAX, MMAX
USE USER_DEFINES
USE NODAL_2D_STORAGE
USE WRITE_DATA
USE POLY_LEVEL_AND_ORDER
USE LOCAL_STORAGE

IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :, :, :) :: EXACT  !< EXACT SOLUTION
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :, :, :) :: ERROR  !< ERRORS AT COLLOCATION POINTS
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: L2_NORM  !< WE EVALUATE THE ERRORS AS L2 NORMS

CONTAINS

!-----------------------------------------------------------------------
! OUTPUT: EXACT, ERROR, L2_NORM
!-----------------------------------------------------------------------
SUBROUTINE GET_ERROR

    IMPLICIT NONE
    
    INTEGER :: I, J, K, S, IERROR
    INTEGER :: PORDERX, PORDERY
    
    DOUBLE PRECISION :: DEL_X, DEL_Y    ! ELEMENT LENGTH
    DOUBLE PRECISION :: L2_NORM_TOTAL(NUM_OF_EQUATION)   ! TOTAL L2 NORM
    
    !-------------------------------------------------------------------
    ALLOCATE(EXACT(0:NMAX, 0:MMAX, NUM_OF_EQUATION, 0:LOCAL_ELEM_NUM-1))
    ALLOCATE(ERROR(0:NMAX, 0:MMAX, NUM_OF_EQUATION, 0:LOCAL_ELEM_NUM-1))
    ALLOCATE(L2_NORM(NUM_OF_EQUATION))
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    EXACT = 0.0D0
    ERROR = 0.0D0
    L2_NORM = 0.0D0
    L2_NORM_TOTAL = 0.0D0
    !-------------------------------------------------------------------
    
    ! GET EXACT SOLUTION------------------------------------------------
    DO K = 0, LOCAL_ELEM_NUM-1
    
        CALL POLY_LEVEL_TO_ORDER(N, PLEVEL_X(K), PORDERX)
        CALL POLY_LEVEL_TO_ORDER(M, PLEVEL_Y(K), PORDERY)
        
        DEL_X = X_LOCAL(4, K) - X_LOCAL(1, K)
        DEL_Y = Y_LOCAL(2, K) - Y_LOCAL(1, K)
        
        CALL EXACT_SOLUTION_GAUSSIAN(PORDERX, PORDERY, &
                                        NUM_OF_EQUATION,&
                                        GL_POINT_X_T(0:PORDERX, PLEVEL_X(K)), &
                                        GL_POINT_Y_T(0:PORDERY, PLEVEL_Y(K)), &
                                        X_LOCAL(1, K), &
                                        Y_LOCAL(1, K), &
                                        DEL_X, DEL_Y, &
                                        EXACT(0:PORDERX, 0:PORDERY, :, K), &
                                        T_TOTAL)
                                                                    
        
        DO S=1, NUM_OF_EQUATION                            
            DO J=0, PORDERY
                DO I=0, PORDERX
                    ERROR(I, J, S, K) = DABS(SOLUTION(I, J, S, K) &
                                        - EXACT(I, J, S, K))
                    L2_NORM(S) = L2_NORM(S) + (ERROR(I, J, S, K))**2
                ENDDO
            
            ENDDO
        ENDDO
        
    ENDDO
    !-------------------------------------------------------------------
    

    CALL MPI_REDUCE(L2_NORM, L2_NORM_TOTAL, NUM_OF_EQUATION, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                    MPI_COMM_WORLD, IERROR)
    
    IF (RANK == 0) THEN
        DO S = 1, NUM_OF_EQUATION
            L2_NORM(S) = DSQRT(L2_NORM_TOTAL(S))
        ENDDO
    ENDIF
    
!    CALL WRITE_ERROR(ERROR)
!    CALL WRITE_RESULTS(N, M, EXACT(:, :, 2), SOLUTION(:, :, 2))
    IF (RANK == 0) THEN
        print *, "L2(1)", L2_NORM(1) 
        print *, "L2(2)", L2_NORM(2) 
        print *, "L2(3)", L2_NORM(3) 
    ENDIF
    
END SUBROUTINE GET_ERROR


END MODULE VERIFICATION
