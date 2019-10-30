!-----------------------------------------------------------------------
!> @brief
!> This module construction the element-based MPI boundaries.
!! Input an element -- if on the MPI boudary -- allocate ghost layer storage
!! -- create MPI WINDOWS -- put solutions on the element interfaces inside
!! the WINDOWS.
!-----------------------------------------------------------------------

MODULE MPI_BOUNDARY

USE MPI
USE hilbert
USE PARAM
USE LOCAL_STORAGE

IMPLICIT NONE 

CONTAINS

!-----------------------------------------------------------------------
!> CONSTRUCT MPI BOUNDARY OF THE LEFT AND RIGHT INTERFACES OF AN ELEMENT
!! (X DIRECTION)
!-----------------------------------------------------------------------
SUBROUTINE MPI_BOUNDARY_CONSTRUCTOR_X(ELEM_K)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: ELEM_K   !< ROOT ELEMENT NUMBER
    INTEGER :: I, J ! DUAL COORD
    INTEGER :: N_NEIGHBOUR  ! NORTH NEIGHTBOUR NUMBER
    INTEGER :: S_NEIGHBOUR  ! SOUTH NEIGHTBOUR NUMBER
    INTEGER :: L_BOUND, R_BOUND ! LOCAL ELEMENT STORAGE BOUNDS
    
    INTEGER(KIND=MPI_ADDRESS_KIND) :: WIN_SIZE, BASEPTR
    INTEGER :: INT_SIZE ! INTEGER SIZE BY BYTES
    INTEGER :: DOUBLE_SIZE ! DOUBLE PRECISION FLOAT NUMBER SIZE BY BYTES
    INTEGER :: IERROR
    
    ! GET DUAL COORD
    CALL d2xy ( EXP_X, ELEM_K, J, I )
    
    ! GET NEIGHBOURS ROOT INDEICES
    CALL xy2d ( EXP_X, J, I-1, S_NEIGHBOUR ) 
    CALL xy2d ( EXP_X, J, I+1, N_NEIGHBOUR ) 
    
    ! LOCAL ELEMENT STORAGE BOUNDARIES
    L_BOUND = ELEM_K 
    R_BOUND = ELEM_K + LOCAL_ELEM_NUM -1
    
    ! Returns the number of bytes occupied by entries in a data type. 
    CALL MPI_TYPE_SIZE(MPI_INTEGER, INT_SIZE, IERROR)
    CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DOUBLE_SIZE, IERROR)
    
    IF ( S_NEIGHBOUR < L_BOUND .OR. S_NEIGHBOUR > R_BOUND) THEN
        ! SOUTH NEIGHTBOUR IS NOT LOCATE LOCALLY, BUILD MPI BOUNDARY
!        CALL MPI_WIN_CREATE(, SIZE, DISP_UNIT, INFO, COMM, WIN, IERROR)
    ENDIF

END SUBROUTINE MPI_BOUNDARY_CONSTRUCTOR_X

SUBROUTINE CREATE_WINDOW(N_MAX)
    
    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: N_MAX    !< MAX POLY ORDER
    
    INTEGER(KIND=MPI_ADDRESS_KIND) :: WIN_SIZE  ! Size of window in bytes (nonnegative integer). 
    INTEGER :: DOUBLE_SIZE ! DOUBLE PRECISION FLOAT NUMBER SIZE BY BYTES
    INTEGER :: WIN
    
    ! Returns the number of bytes occupied by entries in a data type. 
    CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DOUBLE_SIZE, IERROR)
    
    WIN_SIZE = DOUBLE_SIZE * 2 * N_MAX * NUM_OF_EQUATION * LOCAL_ELEM_NUM
    
    CALL MPI_WIN_CREATE(SOLUTION_INT_L, WIN_SIZE, DOUBLE_SIZE, &
                            MPI_INFO_NULL, MPI_COMM_WORLD, &
                            WIN, IERROR)
    
    CALL MPI_WIN_CREATE(SOLUTION_INT_R, WIN_SIZE, DOUBLE_SIZE, &
                            MPI_INFO_NULL, MPI_COMM_WORLD, &
                            WIN, IERROR)
    
END SUBROUTINE CREATE_WINDOW


END MODULE MPI_BOUNDARY
