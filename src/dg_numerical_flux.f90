!-----------------------------------------------------------------------
!> @brief
!> Numerical flux is a combination of two fluxes on the interface.
!! Based on the interface is an element interface or a boundary interface,
!! the numerical flux is computed by either two correspond fluxes on the 
!! element interfaces, or element interface and boundary condition (we 
!! call it external state).
!! On each direction, i.e. x and y, each element has two interfaces. 
!! And the target element's neighbor is decided by its dual coordinate in
!! i-j plane.   
!-----------------------------------------------------------------------

MODULE NUMERICAL_FLUX

USE MPI
USE NODAL_2D_STORAGE
USE INTERFACES_CONSTRUCT
USE PARAM, ONLY: N, M, EXP_X, NUM_OF_EQUATION
USE hilbert
USE RIEMANN_SOLVER
USE AFFINE_MAP
USE POLY_LEVEL_AND_ORDER
USE EXTERNAL_STATE

IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------
! Compute the numerical flux of current element k (direction x)
!-----------------------------------------------------------------------
SUBROUTINE NUMERICAL_FLUX_X(ELEM_K, T)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: ELEM_K   !< ELEMENT NUMBER
    
    INTEGER :: NX   ! POLY ORDER IN X
    INTEGER :: MY   ! POLY ORDER IN Y
    
    INTEGER :: I, J, S
    
    DOUBLE PRECISION :: T   !< CURRENT TIME STEP
    
    DOUBLE PRECISION :: Y
    
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: SOLUTION_EXT ! BOUNDARY CONDITION
    
    ! GET POLY ORDER
    CALL POLY_LEVEL_TO_ORDER(N, PLEVEL_X(ELEM_K), NX)
    CALL POLY_LEVEL_TO_ORDER(M, PLEVEL_Y(ELEM_K), MY)
    
    ! GET DUAL COORDINATE
    CALL d2xy ( EXP_X, ELEM_K, I, J )
    
    !-------------------------------------------------------------------
    IF (I > 0 .AND. I< NUM_OF_ELEMENT_X-1) THEN
    
        ! NOT ON THE BOUNDARY
        CALL RIEMANN1(ELEM_K, I, J, MY)
        
    ELSEIF(I == 0) THEN !-----------------------------------------------
    
        ! ON THE BOTTOM BOUNDARY
        
        ALLOCATE(SOLUTION_EXT(0:MY, NUM_OF_EQUATION))
        SOLUTION_EXT = 0.0D0
        
        DO S = 0, MY
        
            CALL AFFINE_MAPPING(GL_POINT_Y_T(S, PLEVEL_Y(ELEM_K)), &
                                Y, Y_HILBERT(1, ELEM_K), DELTA_Y(ELEM_K))
                                

            CALL EXTERNAL_STATE_GAUSSIAN_EXACT(NUM_OF_EQUATION, &
                                        SOLUTION_EXT(S, :), &
                                         T, X_HILBERT(1, ELEM_K), Y)
                                         
            CALL RIEMANN_X(SOLUTION_EXT(S, :), &
                           SOLUTION_INT_L(S, :, ELEM_K), &
                           NFLUX_X_L(S, :, ELEM_K), -1.0D0)
                           
        ENDDO
    
        DEALLOCATE(SOLUTION_EXT)
        
    ELSEIF(I == NUM_OF_ELEMENT_X-1) THEN !------------------------------
        ! ON THE TOP BOUNDARY

        CALL RIEMANN1(ELEM_K, I, J, MY)
        
        ALLOCATE(SOLUTION_EXT(0:MY, NUM_OF_EQUATION))
        SOLUTION_EXT = 0.0D0
        
        
        DO S = 0, MY
        
            CALL AFFINE_MAPPING(GL_POINT_Y_T(S, PLEVEL_Y(ELEM_K)), &
                                Y, Y_HILBERT(4, ELEM_K), DELTA_Y(ELEM_K))
            
!            CALL EXTERNAL_STATE_GAUSSIAN_EXACT(NUM_OF_EQUATION, &
!                                        SOLUTION_EXT(S, :), &
!                                         T, X_HILBERT(4, ELEM_K), Y)

            ! REFLECT BOUNDARY SOLUTION
            CALL EXTERNAL_STATE_GAUSSIAN_REFLECT(NUM_OF_EQUATION, &
                                            SOLUTION_INT_R(S, :, ELEM_K), &
                                            SOLUTION_EXT(S, :), &
                                            (/1.0D0, 0.0D0/))
                                         
            CALL RIEMANN_X(SOLUTION_INT_R(S, :, ELEM_K), &
                            SOLUTION_EXT(S, :), &
                            NFLUX_X_R(S, :, ELEM_K), 1.0D0)
                            
        ENDDO
        
        DEALLOCATE(SOLUTION_EXT)
    
    ELSE 
        PRINT *, "SOMETHING IS WROING IN 'NUMERICAL_FLUX_X'"
    ENDIF
    !-------------------------------------------------------------------
END SUBROUTINE NUMERICAL_FLUX_X

!-----------------------------------------------------------------------
!> Compute the numerical flux on the left side of the element, and 
!! pass the value of the numerical flux to the corresponding neignbor
!! right boundary with a minus sign. Use this subroutine in x direction.
!-----------------------------------------------------------------------
SUBROUTINE RIEMANN1(ELEM_K, I, J, MY)
    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: ELEM_K
    INTEGER, INTENT(IN) :: I, J
    INTEGER, INTENT(IN) :: MY
    
    INTEGER :: S
    INTEGER :: IDL, IDR ! ELEM NUMBER ON THE TWO SIDE OF THE INTERFACE
    
    IDR = ELEM_K    ! ID ON THE RIGHT SIDE OF THE INTERFACE
        
    CALL xy2d ( EXP_X, I-1, J, IDL )    ! ID ON THE LEFT SIDE OF THE INTERFACE
    
    DO S = 0, MY
        CALL RIEMANN_X(SOLUTION_INT_R(S, :, IDL), &
                       SOLUTION_INT_L(S, :, IDR), &
                       NFLUX_X_L(S, :, IDR), -1.0D0)
                       
        NFLUX_X_R(S, :, IDL) = - NFLUX_X_L(S, :, IDR)
        
    ENDDO
    
END SUBROUTINE RIEMANN1


!-----------------------------------------------------------------------
! Compute the numerical flux of current element k (direction y)
!-----------------------------------------------------------------------
SUBROUTINE NUMERICAL_FLUX_Y(ELEM_K, T)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: ELEM_K   !< ELEMENT NUMBER
    
    INTEGER :: NX   ! POLY ORDER IN X
    INTEGER :: MY   ! POLY ORDER IN Y
    
    INTEGER :: I, J, S
    
    DOUBLE PRECISION :: T   !< CURRENT TIME STEP
    
    DOUBLE PRECISION :: X
    
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: SOLUTION_EXT ! BOUNDARY CONDITION
    
    ! GET POLY ORDER
    CALL POLY_LEVEL_TO_ORDER(N, PLEVEL_X(ELEM_K), NX)
    CALL POLY_LEVEL_TO_ORDER(M, PLEVEL_Y(ELEM_K), MY)
    
    ! GET DUAL COORDINATE
    CALL d2xy ( EXP_X, ELEM_K, I, J )
    
    !-------------------------------------------------------------------
    IF (J > 0 .AND. J< NUM_OF_ELEMENT_Y-1) THEN
    
        ! NOT ON THE BOUNDARY
        CALL RIEMANN2(ELEM_K, I, J, NX)
        
    ELSEIF(J == 0) THEN !-----------------------------------------------
        ! ON THE BOTTOM BOUNDARY
        
        ALLOCATE(SOLUTION_EXT(0:NX, NUM_OF_EQUATION))
        SOLUTION_EXT = 0.0D0
        
        
        DO S = 0, NX
        
            CALL AFFINE_MAPPING(GL_POINT_X_T(S, PLEVEL_X(ELEM_K)), &
                                X, X_HILBERT(1, ELEM_K), DELTA_X(ELEM_K))
        
            CALL EXTERNAL_STATE_GAUSSIAN_EXACT(NUM_OF_EQUATION, &
                                        SOLUTION_EXT(S, :), &
                                         T, X, Y_HILBERT(1, ELEM_K))

            CALL RIEMANN_Y(SOLUTION_EXT(S, :), &
                           SOLUTION_INT_L(S, :, ELEM_K), &
                           NFLUX_Y_D(S, :, ELEM_K), -1.0D0)
        ENDDO
    
        DEALLOCATE(SOLUTION_EXT)
        
    ELSEIF(J == NUM_OF_ELEMENT_Y-1) THEN !------------------------------
    
        ! ON THE TOP BOUNDARY
        
        CALL RIEMANN2(ELEM_K, I, J, NX)
        
        ALLOCATE(SOLUTION_EXT(0:NX, NUM_OF_EQUATION))
        SOLUTION_EXT = 0.0D0
        
        DO S = 0, NX
        
            CALL AFFINE_MAPPING(GL_POINT_X_T(S, PLEVEL_X(ELEM_K)), &
                                X, X_HILBERT(2, ELEM_K), DELTA_X(ELEM_K))
            
            ! EXACT SOLUTION
!            CALL EXTERNAL_STATE_GAUSSIAN_EXACT(NUM_OF_EQUATION, &
!                                        SOLUTION_EXT(S, :), &
!                                         T, X, Y_HILBERT(2, ELEM_K))
            
            ! REFLECT BOUNDARY SOLUTION
            CALL EXTERNAL_STATE_GAUSSIAN_REFLECT(NUM_OF_EQUATION, &
                                            SOLUTION_INT_R(S, :, ELEM_K), &
                                            SOLUTION_EXT(S, :), &
                                            (/0.0D0, 1.0D0/))
                                         
        
            CALL RIEMANN_Y(SOLUTION_INT_R(S, :, ELEM_K), &
                            SOLUTION_EXT(S, :), &
                            NFLUX_Y_U(S, :, ELEM_K), 1.0D0)
                                    
        ENDDO
        
        DEALLOCATE(SOLUTION_EXT)
    
    ELSE 
        PRINT *, "SOMETHING IS WROING IN 'NUMERICAL_FLUX_Y'"
    
    ENDIF
    !-------------------------------------------------------------------
END SUBROUTINE NUMERICAL_FLUX_Y

!-----------------------------------------------------------------------
!> Compute the numerical flux on the left side of the element, and 
!! pass the value of the numerical flux to the corresponding neignbor
!! right boundary with a minus sign. Use this subroutine in y direction.
!-----------------------------------------------------------------------
SUBROUTINE RIEMANN2(ELEM_K, I, J, MX)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: ELEM_K
    INTEGER, INTENT(IN) :: I, J
    INTEGER, INTENT(IN) :: MX
    
    INTEGER :: S
    INTEGER :: IDL, IDR ! ELEM NUMBER ON THE TWO SIDE OF THE INTERFACE
    
    IDR = ELEM_K    ! ID ON THE RIGHT SIDE OF THE INTERFACE
        
    CALL xy2d ( EXP_X, I, J-1, IDL )    ! ID ON THE LEFT SIDE OF THE INTERFACE

    DO S = 0, MX
        CALL RIEMANN_Y(SOLUTION_INT_R(S, :, IDL), &
                       SOLUTION_INT_L(S, :, IDR), &
                       NFLUX_Y_D(S, :, IDR), -1.0D0)
                       
                       
        NFLUX_Y_U(S, :, IDL) = - NFLUX_Y_D(S, :, IDR)
        
    ENDDO
    
END SUBROUTINE RIEMANN2


END MODULE NUMERICAL_FLUX
