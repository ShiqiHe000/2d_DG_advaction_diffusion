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
USE NODAL_2D_STORAGE

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
        
!        if(RANK == 1) THEN
!            PRINT *, ROOT_NUM
!        ENDIF
    
        ! GET DUAL COORD
        CALL d2xy ( EXP_X, ROOT_NUM, J, I )
        
        ! LOCAL ELEMENT STORAGE BOUNDARIES
        L_BOUND = ELEM_RANGE(RANK) + 1
        R_BOUND = ELEM_RANGE(RANK + 1) 
        
        ! X-------------------------------------------------------------
        
        IF(I > 0 .AND. I < NUM_OF_ELEMENT_X-1) THEN     ! NOT ON THE BOUNDARY
        
            CALL CHECK_NORTH(K, I, J, N_NEIGHBOUR, L_BOUND, R_BOUND)
            CALL CHECK_SOUTH(K, I, J, S_NEIGHBOUR, L_BOUND, R_BOUND)

        ELSEIF(I == 0) THEN     ! ON THE BOTTOM, ONLY CHECK NORTH
        
            CALL CHECK_NORTH(K, I, J, N_NEIGHBOUR, L_BOUND, R_BOUND)
        
        ELSEIF(I == NUM_OF_ELEMENT_X - 1) THEN  ! ON THE TOP, ONLY CHECK SOUTH
        
            CALL CHECK_SOUTH(K, I, J, S_NEIGHBOUR, L_BOUND, R_BOUND)
        ELSE 
            PRINT *, "ERROR : SOMETHING WRONG WITH ELEMENT COORDINATES I. SUBROUTINE MPI_BOUNDARY_FLAG"
        ENDIF
                
        !---------------------------------------------------------------
        
        ! Y-------------------------------------------------------------
        
        IF (J > 0 .AND. J < NUM_OF_ELEMENT_Y - 1) THEN
        
            CALL CHECK_WEST(K, I, J, W_NEIGHBOUR, L_BOUND, R_BOUND)
            CALL CHECK_EAST(K, I, J, E_NEIGHBOUR, L_BOUND, R_BOUND)
            
        ELSEIF(J == 0) THEN
        
            CALL CHECK_EAST(K, I, J, E_NEIGHBOUR, L_BOUND, R_BOUND)
        
        ELSEIF(J == NUM_OF_ELEMENT_Y - 1) THEN
        
            CALL CHECK_WEST(K, I, J, W_NEIGHBOUR, L_BOUND, R_BOUND)
            
        ELSE 
            PRINT *, "ERROR : SOMETHING WRONG WITH ELEMENT COORDINATES J. SUBROUTINE MPI_BOUNDARY_FLAG"
        ENDIF 
        !---------------------------------------------------------------
        
    ENDDO
END SUBROUTINE MPI_BOUNDARY_FLAG

SUBROUTINE CHECK_NORTH(K, I, J, NORTH_INDEX, L_BOUND, R_BOUND)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: K   !< ELEMENT NUMBER (LOCAL)
    INTEGER, INTENT(IN) :: I, J     !< ELEMENT COORDINATE
    INTEGER, INTENT(OUT) :: NORTH_INDEX !< NORTH NEIGHBOUR INDEX
    INTEGER, INTENT(IN) :: L_BOUND, R_BOUND !< LOCAL STORAGE BOUNDARY
    
    CALL xy2d ( EXP_X, J, I+1, NORTH_INDEX) 
    
    IF (NORTH_INDEX < L_BOUND .OR. NORTH_INDEX > R_BOUND) THEN
            ! NORTH NEIGHTBOUR IS NOT LOCATE LOCALLY, BUILD MPI BOUNDARY
            MPI_B_FLAG(2, K) = .TRUE.
    ENDIF

END SUBROUTINE CHECK_NORTH

SUBROUTINE CHECK_SOUTH(K, I, J, SOUTH_INDEX, L_BOUND, R_BOUND)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: K   !< ELEMENT NUMBER (LOCAL)
    INTEGER, INTENT(IN) :: I, J     !< ELEMENT COORDINATE
    INTEGER, INTENT(OUT) :: SOUTH_INDEX !< SOUTH NEIGHBOUR INDEX
    INTEGER, INTENT(IN) :: L_BOUND, R_BOUND !< LOCAL STORAGE BOUNDARY
    
    CALL xy2d ( EXP_X, J, I-1, SOUTH_INDEX) 
    
    
    IF ( SOUTH_INDEX < L_BOUND .OR. SOUTH_INDEX > R_BOUND) THEN
        ! SOUTH NEIGHTBOUR IS NOT LOCATE LOCALLY, BUILD MPI BOUNDARY
        MPI_B_FLAG(1, K) = .TRUE.
    ENDIF

END SUBROUTINE CHECK_SOUTH

SUBROUTINE CHECK_WEST(K, I, J, WEST_INDEX, L_BOUND, R_BOUND)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: K   !< ELEMENT NUMBER (LOCAL)
    INTEGER, INTENT(IN) :: I, J     !< ELEMENT COORDINATE
    INTEGER, INTENT(OUT) :: WEST_INDEX !< SOUTH NEIGHBOUR INDEX
    INTEGER, INTENT(IN) :: L_BOUND, R_BOUND !< LOCAL STORAGE BOUNDARY
    
    CALL xy2d ( EXP_X, J-1, I, WEST_INDEX) 
    
     IF ( WEST_INDEX< L_BOUND .OR. WEST_INDEX > R_BOUND) THEN
            ! WEST NEIGHTBOUR IS NOT LOCATE LOCALLY, BUILD MPI BOUNDARY
            MPI_B_FLAG(3, K) = .TRUE.
    ENDIF

END SUBROUTINE CHECK_WEST

SUBROUTINE CHECK_EAST(K, I, J, EAST_INDEX, L_BOUND, R_BOUND)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: K   !< ELEMENT NUMBER (LOCAL)
    INTEGER, INTENT(IN) :: I, J     !< ELEMENT COORDINATE
    INTEGER, INTENT(OUT) :: EAST_INDEX !< SOUTH NEIGHBOUR INDEX
    INTEGER, INTENT(IN) :: L_BOUND, R_BOUND !< LOCAL STORAGE BOUNDARY
    
    CALL xy2d ( EXP_X, J+1, I, EAST_INDEX) 
    
     IF (EAST_INDEX < L_BOUND .OR. EAST_INDEX > R_BOUND) THEN
        ! EAST NEIGHTBOUR IS NOT LOCATE LOCALLY, BUILD MPI BOUNDARY
        MPI_B_FLAG(4, K) = .TRUE.
    ENDIF

END SUBROUTINE CHECK_EAST

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
