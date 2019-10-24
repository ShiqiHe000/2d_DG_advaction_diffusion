!-----------------------------------------------------------------------
!> @brief
!> In order to get time ferivative, we need spatial derivative first.
!! Here we compute coefficient matrix A times spactial derivative together.
!-----------------------------------------------------------------------

MODULE A_TIMES_SPATIAL_DER

USE MPI
USE POLY_LEVEL_AND_ORDER
USE PARAM, ONLY: N, M, NUM_OF_EQUATION
USE NODAL_2D_STORAGE
USE FLUX_VECTORS
USE SPACTIAL_DERIVATIVE
USE LOCAL_STORAGE

IMPLICIT NONE

CONTAINS

SUBROUTINE A_TIMES_SPATIAL_DERIRATIVE_X(ELEM_K)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: ELEM_K   !< CURRENT ELEMENT K
    INTEGER :: PORDERX, PORDERY   ! POLY ORDER
    INTEGER :: I, J, S
    
    CALL POLY_LEVEL_TO_ORDER(N, PLEVEL_X(ELEM_K), PORDERX)
    CALL POLY_LEVEL_TO_ORDER(M, PLEVEL_Y(ELEM_K), PORDERY)
    
    ! A * F ------------------------------------------------------------
    DO J = 0, PORDERY
        DO I = 0, PORDERX
            CALL XFLUX(SOLUTION(I, J, :, ELEM_K), FLUX_X(I, J, :, ELEM_K))
        ENDDO
        
        ! A * F'--------------------------------------------------------
        CALL DG_SPATIAL_DERIVATIVE(PORDERX, NFLUX_X_L(J, :, ELEM_K), &
                                    NFLUX_X_R(J, :, ELEM_K), &
                                    FLUX_X(0:PORDERX, J, :, ELEM_K), &
                                    FLUX_DER_X(0:PORDERX, J, :, ELEM_K), &
                                    M_FIRST_DER_X_T(0:PORDERX, 0:PORDERX, PLEVEL_X(ELEM_K)), &
                                    LAGRANGE_LEFT_T(0:PORDERX, PLEVEL_X(ELEM_K)), &
                                    LAGRANGE_RIGHT_T(0:PORDERX, PLEVEL_X(ELEM_K)), &
                                    GL_W_X_T(0:PORDERX, PLEVEL_X(ELEM_K)))

        DO S = 1, NUM_OF_EQUATION
            DO I = 0, PORDERX
                SOLUTION_TIME_DER(I, J, S, ELEM_K) = &
                                                - (2.0D0 / DELTA_X(ELEM_K)) &
                                                * FLUX_DER_X(I, J, S, ELEM_K)
            ENDDO
        
        ENDDO
        !---------------------------------------------------------------
    ENDDO
    !-------------------------------------------------------------------

END SUBROUTINE A_TIMES_SPATIAL_DERIRATIVE_X

SUBROUTINE A_TIMES_SPATIAL_DERIRATIVE_Y(ELEM_K)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: ELEM_K   !< CURRENT ELEMENT K
    INTEGER :: PORDERX, PORDERY   ! POLY ORDER
    INTEGER :: I, J, S
    
    CALL POLY_LEVEL_TO_ORDER(N, PLEVEL_X(ELEM_K), PORDERX)
    CALL POLY_LEVEL_TO_ORDER(M, PLEVEL_Y(ELEM_K), PORDERY)
    
    ! B * F ------------------------------------------------------------
    DO I = 0, PORDERX
        DO J = 0, PORDERY
            CALL YFLUX(SOLUTION(I, J, :, ELEM_K), FLUX_Y(I, J, :, ELEM_K))
        ENDDO
        
        ! B * F'--------------------------------------------------------
        CALL DG_SPATIAL_DERIVATIVE(PORDERY, NFLUX_Y_D(I, :, ELEM_K), &
                                    NFLUX_Y_U(I, :, ELEM_K), &
                                    FLUX_Y(I, 0:PORDERY, :, ELEM_K), &
                                    FLUX_DER_Y(I, 0:PORDERY, :, ELEM_K), &
                                    M_FIRST_DER_Y_T(0:PORDERY, 0:PORDERY, PLEVEL_Y(ELEM_K)), &
                                    LAGRANGE_DOWN_T(0:PORDERY, PLEVEL_Y(ELEM_K)), &
                                    LAGRANGE_UP_T(0:PORDERY, PLEVEL_Y(ELEM_K)), &
                                    GL_W_Y_T(0:PORDERY, PLEVEL_Y(ELEM_K)))
        
        DO S = 1, NUM_OF_EQUATION
            DO J = 0, PORDERY
                SOLUTION_TIME_DER(I, J, S, ELEM_K) = SOLUTION_TIME_DER(I, J, S, ELEM_K) &
                                                - (2.0D0 / DELTA_Y(ELEM_K)) &
                                                * FLUX_DER_Y(I, J, S, ELEM_K)
            ENDDO
        
        ENDDO
        !---------------------------------------------------------------
    ENDDO
    !-------------------------------------------------------------------
    

END SUBROUTINE A_TIMES_SPATIAL_DERIRATIVE_Y

END MODULE A_TIMES_SPATIAL_DER
