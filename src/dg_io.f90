!-----------------------------------------------------------------------
!> @brief
!> This module contains serial IO and parallel IO
!-----------------------------------------------------------------------

MODULE IO

USE MPI
USE LOCAL_STORAGE
USE PARAM
USE OUTPUT

IMPLICIT NONE

CONTAINS

SUBROUTINE SERIAL_IO(T)

    IMPLICIT NONE 
    
!    INTEGER, INTENT(IN) :: TOTAL_ELEM   !< TOTAL ELEMENT NUMBER, START WITH 1
!    INTEGER :: RECVCOUNTS(NUM_PROC)     ! containing the number of elements that are received from each process 
!    INTEGER :: DISPLS1(NUM_PROC) ! Entry i specifies the displacement relative to recvbuf at which to place the incoming data from process i 
    
    INTEGER :: I
    
    DOUBLE PRECISION, INTENT(IN) :: T   !< CURRENT TIME 
    
!    DOUBLE PRECISION :: X_GLOBAL(4, TOTAL_ELEM) ! GLOBAL ELEM COORDS
!    DOUBLE PRECISION :: Y_GLOBAL(4, TOTAL_ELEM)
!    DOUBLE PRECISION :: SOLUTION_ALL(0:NMAX, 0:MMAX, NUM_OF_EQUATION, 0:TOTAL_ELEM-1)
    
!    RECVCOUNTS = 0
!    DISPLS = 0
!    X_GLOBAL = 0.0D0, Y_GLOBAL = 0.0D0
!    SOLUTION_ALL = 0.0D0
    
!    IF (RANK == 0) THEN
!        RECVCOUNTS(1) = ELEM_RANGE(1) + 1
    
!        DO I = 2, NUM_PROC
    
!            RECVCOUNTS(I) = ELEM_RANGE(I) - ELEM_RANGE(I-1)
!            DISPLS1(I) = DISPLS(I-1) + RECVCOUNTS(I-1) * 4
            
!        ENDDO
!    ENDIF
    
    ! GATHER INFORMATION------------------------------------------------
!    CALL MPI_GATHERV(X_LOCAL, LOCAL_ELEM_NUM*4, MPI_DOUBLE_PRECISION, &
!                    X_GLOBAL, RECVCOUNTS*4, DISPLS1, MPI_DOUBLE_PRECISION, &
!                    & 0, MPI_COMM_WORLD, IERROR)
                    
!    CALL MPI_GATHERV(Y_LOCAL, LOCAL_ELEM_NUM*4, MPI_DOUBLE_PRECISION, &
!                    Y_GLOBAL, RECVCOUNTS*4, DISPLS1, MPI_DOUBLE_PRECISION, &
!                    & 0, MPI_COMM_WORLD, IERROR)
    
!    CALL MPI_GATHERV()
    !-------------------------------------------------------------------
    
!    MPI_GATHERV(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNTS,
!        DISPLS, RECVTYPE, ROOT, COMM, IERROR)
    
    ! ONLY PROC OUTPUT RESULT SEQUENTIALLY
    DO I = 0, NUM_PROC-1
        IF(RANK == I) THEN
            CALL WRITE_MESH(LOCAL_ELEM_NUM, X_LOCAL, Y_LOCAL, &
                            PLEVEL_X, PLEVEL_Y, &
                            SOLUTION, T)
        
        ELSE
            CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
        
        ENDIF
    ENDDO


END SUBROUTINE SERIAL_IO


END MODULE IO
