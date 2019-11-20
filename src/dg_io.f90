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
    
    INTEGER :: I
    
    DOUBLE PRECISION, INTENT(IN) :: T   !< CURRENT TIME 
    
    ! ONLY PROC OUTPUT RESULT SEQUENTIALLY
    DO I = 0, NUM_PROC-1
!    DO I = 0, 2
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
