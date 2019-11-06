!-----------------------------------------------------------------------
!> @brief
!> This module construction the element-based MPI boundaries.
!! Flag boundary elements -- create windows for fluxes exchange
!-----------------------------------------------------------------------

MODULE MPI_BOUNDARY

USE MPI
USE hilbert
USE PARAM
USE LOCAL_STORAGE

IMPLICIT NONE 

CONTAINS

!-----------------------------------------------------------------------
!> FLAG THE ELEMENT WHO HAS A MPI BOUNDARY.
!! ASSIGN THE FLAG = .TRUE. TO THE CORRSPONDING ELEMENT FACES.
!! THE ELEMENT FACES ARE DEFINES BELOW
!!                f2
!!          -------------
!!          |           |
!!          |           |
!!     f3   |           |   f4
!!          |           |
!!          -------------
!!                f1
!-----------------------------------------------------------------------
SUBROUTINE MPI_BOUNDARY_FLAG

    IMPLICIT NONE 
    
    INTEGER :: K
    INTEGER :: I, J ! DUAL COORD
    INTEGER :: N_NEIGHBOUR  ! NORTH NEIGHTBOUR NUMBER
    INTEGER :: S_NEIGHBOUR  ! SOUTH NEIGHTBOUR NUMBER
    INTEGER :: W_NEIGHBOUR  ! WEST NEIGHTBOUR NUMBER
    INTEGER :: E_NEIGHBOUR  ! EAST NEIGHTBOUR NUMBER
    INTEGER :: L_BOUND, R_BOUND ! LOCAL ELEMENT STORAGE BOUNDS
    INTEGER :: ROOT_NUM ! ROOT ELEMENT NUMBER
    
    DO K = 0, LOCAL_ELEM_NUM-1
    
        ROOT_NUM = K + ELEM_RANGE(RANK) + 1
    
        ! GET DUAL COORD
        CALL d2xy ( EXP_X, ROOT_NUM, J, I )
        
        ! LOCAL ELEMENT STORAGE BOUNDARIES
        L_BOUND = ELEM_RANGE(RANK) + 1
        R_BOUND = ELEM_RANGE(RANK + 1) 
        
        ! X-------------------------------------------------------------
        
        ! GET NEIGHBOURS ROOT INDEICES
        CALL xy2d ( EXP_X, J, I-1, S_NEIGHBOUR ) 
        CALL xy2d ( EXP_X, J, I+1, N_NEIGHBOUR ) 
        
        IF ( S_NEIGHBOUR < L_BOUND .OR. S_NEIGHBOUR > R_BOUND) THEN
            ! SOUTH NEIGHTBOUR IS NOT LOCATE LOCALLY, BUILD MPI BOUNDARY
            MPI_B_FLAG(1, K) = .TRUE.
        ENDIF
        
        IF (N_NEIGHBOUR < L_BOUND .OR. N_NEIGHBOUR > R_BOUND) THEN
            ! NORTH NEIGHTBOUR IS NOT LOCATE LOCALLY, BUILD MPI BOUNDARY
            MPI_B_FLAG(2, K) = .TRUE.
        ENDIF
            
        !---------------------------------------------------------------
        
        ! Y-------------------------------------------------------------
        
        ! GET NEIGHBOURS ROOT INDEICES
        CALL xy2d ( EXP_X, J-1, I, W_NEIGHBOUR ) 
        CALL xy2d ( EXP_X, J+1, I, E_NEIGHBOUR ) 
        
        IF ( W_NEIGHBOUR < L_BOUND .OR. W_NEIGHBOUR > R_BOUND) THEN
            ! WEST NEIGHTBOUR IS NOT LOCATE LOCALLY, BUILD MPI BOUNDARY
            MPI_B_FLAG(3, K) = .TRUE.
        ENDIF
        
        IF (E_NEIGHBOUR < L_BOUND .OR. E_NEIGHBOUR > R_BOUND) THEN
            ! EAST NEIGHTBOUR IS NOT LOCATE LOCALLY, BUILD MPI BOUNDARY
            MPI_B_FLAG(4, K) = .TRUE.
        ENDIF
            
        !---------------------------------------------------------------
        
    ENDDO
END SUBROUTINE MPI_BOUNDARY_FLAG


!-----------------------------------------------------------------------
!> Create windows to have remotely accessible
!! memories.
!-----------------------------------------------------------------------
SUBROUTINE CREATE_WINDOW(COLUMN, WIN1, WIN2, Q1, Q2)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: COLUMN   !< INPUT VARIABLE FIRST COLUMN SIZE, START FROM 0
    
    INTEGER, INTENT(OUT) :: WIN1, WIN2  !< Window object returned by the call (handle)
    INTEGER(KIND=MPI_ADDRESS_KIND) :: WIN_SIZE
    INTEGER :: DOUBLE_SIZE ! DOUBLE PRECISION FLOAT NUMBER SIZE BY BYTES
    INTEGER :: DISP_UNIT    ! Local unit size for displacements, in bytes (positive integer)
    INTEGER :: IERROR
    
    DOUBLE PRECISION :: Q1(0:COLUMN, NUM_OF_EQUATION, 0:LOCAL_ELEM_NUM-1)    !< Initial address of window1
    DOUBLE PRECISION :: Q2(0:COLUMN, NUM_OF_EQUATION, 0:LOCAL_ELEM_NUM-1)    !< Initial address of window2
    
    ! Returns the number of bytes occupied by entries in a data type. 
    CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DOUBLE_SIZE, IERROR)
    
    DISP_UNIT = DOUBLE_SIZE
    
    WIN_SIZE = DOUBLE_SIZE * (COLUMN + 1) * NUM_OF_EQUATION * LOCAL_ELEM_NUM
    
    CALL MPI_WIN_CREATE(Q1, WIN_SIZE, DISP_UNIT, MPI_INFO_NULL, &
                        MPI_COMM_WORLD, WIN1, IERROR)
    CALL MPI_WIN_CREATE(Q2, WIN_SIZE, DISP_UNIT, MPI_INFO_NULL, &
                        MPI_COMM_WORLD, WIN2, IERROR)

END SUBROUTINE CREATE_WINDOW


END MODULE MPI_BOUNDARY