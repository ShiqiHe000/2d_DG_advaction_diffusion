!-----------------------------------------------------------------------
!> @brief
!> Numerical flux is a combination of two fluxes on the interface.
!! Based on the interface is either an element interface or a boundary interface,
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
USE LOCAL_STORAGE
USE SEARCH_RANK
USE INDEX_LOCAL_GLOBAL

IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------
! Compute the numerical flux of current element k (direction x)
!-----------------------------------------------------------------------
SUBROUTINE NUMERICAL_FLUX_X(LELEM_K, T)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: LELEM_K   !< ELEMENT NUMBER (LOCAL)
    INTEGER :: RELEM_K   ! ELEMENT NUMBER (ROOT)
    
    INTEGER :: NX   ! POLY ORDER IN X
    INTEGER :: MY   ! POLY ORDER IN Y
    
    INTEGER :: I, J, S
    
    DOUBLE PRECISION :: DEL_Y    ! ELEMENT LENGTH IN Y DIRECTION
    
    DOUBLE PRECISION :: T   !< CURRENT TIME STEP
    
    DOUBLE PRECISION :: Y
    
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: SOLUTION_EXT ! BOUNDARY CONDITION
    
    ! GET POLY ORDER
    CALL POLY_LEVEL_TO_ORDER(N, PLEVEL_X(LELEM_K), NX)
    CALL POLY_LEVEL_TO_ORDER(M, PLEVEL_Y(LELEM_K), MY)
    
    ! GET ELEMENT ROOT NUMBER
    CALL INDEX_LOCAL_TO_GLOBAL(RANK, LELEM_K, RELEM_K)
    
    ! GET DUAL COORDINATE
    CALL d2xy ( EXP_X, RELEM_K, J, I )
    
    !-------------------------------------------------------------------
    IF (I > 0 .AND. I< NUM_OF_ELEMENT_X-1) THEN
!        print *, "I", I, RELEM_K
    
        ! NOT ON THE BOUNDARY
        CALL RIEMANN1(RELEM_K, I, J, MY)
        
    ELSEIF(I == 0) THEN !-----------------------------------------------
    
        ! ON THE BOTTOM BOUNDARY
        
        ALLOCATE(SOLUTION_EXT(0:MY, NUM_OF_EQUATION))
        SOLUTION_EXT = 0.0D0
        
    
        DO S = 0, MY
        
            DEL_Y = Y_LOCAL(2, LELEM_K) - Y_LOCAL(1, LELEM_K)
        
            CALL AFFINE_MAPPING(GL_POINT_Y_T(S, PLEVEL_Y(LELEM_K)), &
                                Y, Y_LOCAL(1, LELEM_K), DEL_Y)
                                

            CALL EXTERNAL_STATE_GAUSSIAN_EXACT(NUM_OF_EQUATION, &
                                        SOLUTION_EXT(S, :), &
                                         T, X_LOCAL(1, LELEM_K), Y)
                                         
            CALL RIEMANN_X(SOLUTION_EXT(S, :), &
                           SOLUTION_INT_L(S, :, LELEM_K), &
                           NFLUX_X_L(S, :, LELEM_K), -1.0D0)
                           
        ENDDO
    
        DEALLOCATE(SOLUTION_EXT)
        
    ELSEIF(I == NUM_OF_ELEMENT_X-1) THEN !------------------------------
        ! ON THE TOP BOUNDARY
        
        ! LEFT INTERFACE
        CALL RIEMANN1(RELEM_K, I, J, MY)
        
        ALLOCATE(SOLUTION_EXT(0:MY, NUM_OF_EQUATION))
        SOLUTION_EXT = 0.0D0
        
        ! RIGHT INTERFACE
        DO S = 0, MY
        
            DEL_Y = Y_LOCAL(3, LELEM_K) - Y_LOCAL(4, LELEM_K)
            
            CALL AFFINE_MAPPING(GL_POINT_Y_T(S, PLEVEL_Y(LELEM_K)), &
                                Y, Y_LOCAL(4, LELEM_K), DEL_Y)
            
            CALL EXTERNAL_STATE_GAUSSIAN_EXACT(NUM_OF_EQUATION, &
                                        SOLUTION_EXT(S, :), &
                                         T, X_LOCAL(4, LELEM_K), Y)

            ! REFLECT BOUNDARY SOLUTION
!            CALL EXTERNAL_STATE_GAUSSIAN_REFLECT(NUM_OF_EQUATION, &
!                                            SOLUTION_INT_R(S, :, LELEM_K), &
!                                            SOLUTION_EXT(S, :), &
!                                            (/1.0D0, 0.0D0/))
                                         
            CALL RIEMANN_X(SOLUTION_INT_R(S, :, LELEM_K), &
                            SOLUTION_EXT(S, :), &
                            NFLUX_X_R(S, :, LELEM_K), 1.0D0)
                            
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
SUBROUTINE RIEMANN1(LELEM_K, I, J, MY)
    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: LELEM_K   !< LOCAL ELEMENT NUMBER
!    INTEGER, INTENT(IN) :: RELEM_K   !< ROOT ELEMENT NUMBER
    INTEGER, INTENT(IN) :: I, J     !< ELEMENT COORDINATE
    INTEGER, INTENT(IN) :: MY       !< ELEMENT LEFT Y INTERFACE POLYNOMIAL ORDER
    
    INTEGER :: S
    INTEGER :: ENTRY_COUNT  ! NUMBER OF DATA TO GET FROM THE REMOTE RANK
    INTEGER :: IDL, IDR ! ELEM NUMBER ON THE TWO SIDE OF THE INTERFACE (LOCAL)
    INTEGER :: IDL_G    ! ELEMENT ON THE LEFT SIDE OF THE INTERFACE: GLOBAL INDEX
    
    INTEGER :: TARGET_RANK  ! THE NEIBOURHOUR ELEMENT'S RANK
    INTEGER(KIND=MPI_ADDRESS_KIND) :: TARGET_DISP   ! Displacement from window start to the beginning of the target buffer
    INTEGER :: ORIGIN_COUNT ! NUMBER OF DATATO PUT ON THE REMOTE RANK
    
    DOUBLE PRECISION :: REMOTE_SOLUTION_INT_L(0:MMAX, NUM_OF_EQUATION)  ! SOLUTION_INT_L FROM REMOTE BUFFER
    
    REMOTE_SOLUTION_INT_L = 0.0D0
    
    IDR = LELEM_K    ! ID ON THE RIGHT SIDE OF THE INTERFACE (LOCAL)
        
    CALL xy2d ( EXP_X, J, I-1, IDL_G )    ! ID ON THE LEFT SIDE OF THE INTERFACE
    
    ! NEED REMOTE INFORMATION
    IF (MPI_B_FLAG(1, LELEM_K)) THEN
    
        ENTRY_COUNT = (MMAX + 1) * NUM_OF_EQUATION
        
        CALL FIND_RANK(IDL_G, TARGET_RANK)
        
        CALL INDEX_GLOBAL_TO_LOCAL(TARGET_RANK, IDL_G, IDL)
        
        TARGET_DISP = IDL * (1 + MMAX) * NUM_OF_EQUATION
        
        CALL MPI_GET(REMOTE_SOLUTION_INT_L, ENTRY_COUNT, &
                        MPI_DOUBLE_PRECISION, TARGET_RANK, &
                        TARGET_DISP, ENTRY_COUNT, MPI_DOUBLE_PRECISION, &
                        WIN_INTERFACE_R, IERROR)
        
        ! NOW WE ASSUMING CONFORMING INTERFACES
        DO S = 0, MY
            CALL RIEMANN_X(REMOTE_SOLUTION_INT_L(S, :), &
                           SOLUTION_INT_L(S, :, IDR), &
                           NFLUX_X_L(S, :, IDR), -1.0D0)
        ENDDO
        
        ORIGIN_COUNT = (MY + 1) * NUM_OF_EQUATION
        
!        TARGET_DISP = 0
               
        CALL MPI_PUT( - NFLUX_X_L(0:MY, :, IDR), ORIGIN_COUNT, &
                    MPI_DOUBLE_PRECISION, TARGET_RANK, &
                    TARGET_DISP, ORIGIN_COUNT, &
                    MPI_DOUBLE_PRECISION, WIN_INTERFACE_L, IERROR)
                
    ELSE ! NEED LOCAL INFORMATION
        CALL INDEX_GLOBAL_TO_LOCAL(RANK, IDL_G, IDL)
    
        DO S = 0, MY
            CALL RIEMANN_X(SOLUTION_INT_R(S, :, IDL), &
                           SOLUTION_INT_L(S, :, IDR), &
                           NFLUX_X_L(S, :, IDR), -1.0D0)
                           
            NFLUX_X_R(S, :, IDL) = - NFLUX_X_L(S, :, IDR)
            
        ENDDO
    ENDIF
    
END SUBROUTINE RIEMANN1


!-----------------------------------------------------------------------
! Compute the numerical flux of current element k (direction y)
!-----------------------------------------------------------------------
SUBROUTINE NUMERICAL_FLUX_Y(LELEM_K, T)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: LELEM_K   !< ELEMENT NUMBER (LOCAL)
    INTEGER :: RELEM_K   !< ELEMENT NUMBER (ROCAL)
    
    INTEGER :: NX   ! POLY ORDER IN X
    INTEGER :: MY   ! POLY ORDER IN Y
    
    INTEGER :: I, J, S
    
    DOUBLE PRECISION :: T   !< CURRENT TIME STEP
    
    DOUBLE PRECISION :: X
    DOUBLE PRECISION :: DEL_X   ! ELEMENT LENGTH IN X DIRECTION
    
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: SOLUTION_EXT ! BOUNDARY CONDITION
    
    ! GET POLY ORDER
    CALL POLY_LEVEL_TO_ORDER(N, PLEVEL_X(LELEM_K), NX)
    CALL POLY_LEVEL_TO_ORDER(M, PLEVEL_Y(LELEM_K), MY)
    
    RELEM_K = LELEM_K + ELEM_RANGE(RANK) + 1
    
    ! GET DUAL COORDINATE
    CALL d2xy ( EXP_X, RELEM_K, J, I )
    
    !-------------------------------------------------------------------
    IF (J > 0 .AND. J< NUM_OF_ELEMENT_Y-1) THEN
    
        ! NOT ON THE BOUNDARY
        CALL RIEMANN2(RELEM_K, I, J, NX)
        
    ELSEIF(J == 0) THEN !-----------------------------------------------
        ! ON THE BOTTOM BOUNDARY
        
        ALLOCATE(SOLUTION_EXT(0:NX, NUM_OF_EQUATION))
        SOLUTION_EXT = 0.0D0
        
        
        DO S = 0, NX
        
            DEL_X = X_LOCAL(3, LELEM_K) - X_LOCAL(1, LELEM_K)
        
            CALL AFFINE_MAPPING(GL_POINT_X_T(S, PLEVEL_X(LELEM_K)), &
                                X, X_LOCAL(1, LELEM_K), DEL_X)
        
            CALL EXTERNAL_STATE_GAUSSIAN_EXACT(NUM_OF_EQUATION, &
                                        SOLUTION_EXT(S, :), &
                                         T, X, Y_LOCAL(1, LELEM_K))

            CALL RIEMANN_Y(SOLUTION_EXT(S, :), &
                           SOLUTION_INT_L(S, :, LELEM_K), &
                           NFLUX_Y_D(S, :, LELEM_K), -1.0D0)
        ENDDO
    
        DEALLOCATE(SOLUTION_EXT)
        
    ELSEIF(J == NUM_OF_ELEMENT_Y-1) THEN !------------------------------
    
        ! ON THE TOP BOUNDARY
        
        CALL RIEMANN2(RELEM_K, I, J, NX)
        
        ALLOCATE(SOLUTION_EXT(0:NX, NUM_OF_EQUATION))
        SOLUTION_EXT = 0.0D0
        
        DO S = 0, NX
        
            DEL_X = X_LOCAL(3, LELEM_K) - X_LOCAL(1, LELEM_K)
        
            CALL AFFINE_MAPPING(GL_POINT_X_T(S, PLEVEL_X(LELEM_K)), &
                                X, X_LOCAL(2, LELEM_K), DEL_X)
            
            ! EXACT SOLUTION
            CALL EXTERNAL_STATE_GAUSSIAN_EXACT(NUM_OF_EQUATION, &
                                        SOLUTION_EXT(S, :), &
                                         T, X, Y_LOCAL(2, LELEM_K))
            
            ! REFLECT BOUNDARY SOLUTION
!            CALL EXTERNAL_STATE_GAUSSIAN_REFLECT(NUM_OF_EQUATION, &
!                                            SOLUTION_INT_R(S, :, LELEM_K), &
!                                            SOLUTION_EXT(S, :), &
!                                            (/0.0D0, 1.0D0/))
                                         
        
            CALL RIEMANN_Y(SOLUTION_INT_R(S, :, LELEM_K), &
                            SOLUTION_EXT(S, :), &
                            NFLUX_Y_U(S, :, LELEM_K), 1.0D0)
                                    
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
    
    INTEGER, INTENT(IN) :: ELEM_K   !< ROOT ELEMENT NUMBER
    INTEGER, INTENT(IN) :: I, J     !< ELEMENT COORDINATE
    INTEGER, INTENT(IN) :: MX       !< ELEMENT LEFT X INTERFACE POLYNOMIAL ORDER
    
    INTEGER :: S, P, TOTAL_CELLS, ENTRY_COUNT
    INTEGER :: IDL, IDR ! ELEM NUMBER ON THE TWO SIDE OF THE INTERFACE
    
    INTEGER :: TARGET_RANK  ! THE NEIBOURHOUR ELEMENT'S RANK
    INTEGER(KIND=MPI_ADDRESS_KIND) :: TARGET_DISP   ! Displacement from window start to the beginning of the target buffer
    INTEGER :: ORIGIN_COUNT ! Number of entries in origin buffer
    
    IDR = ELEM_K    ! ID ON THE RIGHT SIDE OF THE INTERFACE
        
    CALL xy2d ( EXP_X, J-1, I, IDL )    ! ID ON THE LEFT SIDE OF THE INTERFACE

    IF (MPI_B_FLAG(3, ELEM_K)) THEN
    
        P = NMAX+1
        TOTAL_CELLS = (NMAX + 1)*2 - 1
        ENTRY_COUNT = (TOTAL_CELLS - P + 1) * NUM_OF_EQUATION
        TARGET_DISP = P 
        
        CALL FIND_RANK(IDL, TARGET_RANK)
        
        print *, IDL, TARGET_RANK
        
        CALL MPI_GET(SOLUTION_INT_L(P:TOTAL_CELLS, :, IDR), ENTRY_COUNT, &
                        MPI_DOUBLE_PRECISION, TARGET_RANK, &
                        TARGET_DISP, ENTRY_COUNT, MPI_DOUBLE_PRECISION, &
                        WIN_INTERFACE_R, IERROR)

        DO S = 0, MX
            CALL RIEMANN_Y(SOLUTION_INT_R(S, :, IDR), &
                           SOLUTION_INT_L(S, :, IDR), &
                           NFLUX_Y_D(S, :, IDR), -1.0D0)
        ENDDO
        
        ORIGIN_COUNT = (MX + 1) * NUM_OF_EQUATION
        
        TARGET_DISP = 0
               
        CALL MPI_PUT( - NFLUX_Y_D(0:MX, :, IDR), ORIGIN_COUNT, &
                    MPI_DOUBLE_PRECISION, TARGET_RANK, &
                    TARGET_DISP, ORIGIN_COUNT, &
                    MPI_DOUBLE_PRECISION, WIN_INTERFACE_L, IERROR)
                
    ELSE ! NEED LOCAL INFORMATION
        DO S = 0, MX
            CALL RIEMANN_Y(SOLUTION_INT_R(S, :, IDL), &
                           SOLUTION_INT_L(S, :, IDR), &
                           NFLUX_Y_D(S, :, IDR), -1.0D0)
                           
                           
            NFLUX_Y_U(S, :, IDL) = - NFLUX_Y_D(S, :, IDR)
        
        ENDDO
    ENDIF
    
END SUBROUTINE RIEMANN2


END MODULE NUMERICAL_FLUX
