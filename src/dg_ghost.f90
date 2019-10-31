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

!-----------------------------------------------------------------------
!> If element interfaces (x direction) are on the MPI boundary, 
!! then copy the solutions to the ghost layer. 
!-----------------------------------------------------------------------
SUBROUTINE GHOST_CONSTUCT_X(ELEM_K, PORDER_Y)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: ELEM_K   !< LOCAL ELEMENT(ROOT) NUMBER, START WITH 0
    INTEGER, INTENT(IN) :: PORDER_Y   !< POLY ORDER
    INTEGER :: I, S, L
    
    S = MMAX+1
    
    ! FACE 1
    IF(MPI_B_FLAG(1, ELEM_K)) THEN
        DO L = 1, NUM_OF_EQUATION
            DO I = 0, PORDER_Y
                SOLUTION_INT_L(S + I, L, ELEM_K) = SOLUTION_INT_L(I, L, ELEM_K)
            ENDDO
        ENDDO
    ENDIF
    
    ! FACE 2
    IF(MPI_B_FLAG(2, ELEM_K)) THEN
        DO L = 1, NUM_OF_EQUATION
            DO I = 0, PORDER_Y
                SOLUTION_INT_R(S + I, L, ELEM_K) = SOLUTION_INT_R(I, L, ELEM_K)
            ENDDO
        ENDDO
    ENDIF

END SUBROUTINE GHOST_CONSTUCT_X



!-----------------------------------------------------------------------
!> If element interfaces (y direction) are on the MPI boundary, 
!! then copy the solutions to the ghost layer. 
!-----------------------------------------------------------------------
SUBROUTINE GHOST_CONSTUCT_Y(ELEM_K, PORDER_X)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: ELEM_K   !< LOCAL ELEMENT(ROOT) NUMBER, START WITH 0
    INTEGER, INTENT(IN) :: PORDER_X   !< POLY ORDER
    INTEGER :: I, S, L
    
    S = NMAX+1
    
    ! FACE 1
    IF(MPI_B_FLAG(3, ELEM_K)) THEN
        DO L = 1, NUM_OF_EQUATION
            DO I = 0, PORDER_X
                SOLUTION_INT_L(S + I, L, ELEM_K) = SOLUTION_INT_L(I, L, ELEM_K)
            ENDDO
        ENDDO
    ENDIF
    
    ! FACE 2
    IF(MPI_B_FLAG(4, ELEM_K)) THEN
        DO L = 1, NUM_OF_EQUATION
            DO I = 0, PORDER_X
                SOLUTION_INT_R(S + I, L, ELEM_K) = SOLUTION_INT_R(I, L, ELEM_K)
            ENDDO
        ENDDO
    ENDIF

END SUBROUTINE GHOST_CONSTUCT_Y

END MODULE GHOST
