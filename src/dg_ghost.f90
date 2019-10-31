!-----------------------------------------------------------------------
!> @brief
!> Manage ghost cells message exchange
!-----------------------------------------------------------------------

MODULE GHOST

USE MPI
USE LOCAL_STORAGE
USE PARAM

IMPLICIT NONE 

CONTAINS

SUBROUTINE GHOST_CONSTUCT_X(ELEM_K, PORDER_Y)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: ELEM_K   !< LOCAL ELEMENT(ROOT) NUMBER, START WITH 0
    INTEGER, INTENT(IN) :: PORDER_Y   !< POLY ORDER
    INTEGER :: I
    
    ! FACE 1
    IF(MPI_B_FLAG(1, ELEM_K)) THEN
        DO I = 0, PORDER_Y
            SOLUTION_INT_L()
        ENDDO
    ENDIF

END SUBROUTINE GHOST_CONSTUCT_X

END MODULE GHOST
