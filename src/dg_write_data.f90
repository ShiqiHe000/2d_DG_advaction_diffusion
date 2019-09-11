!-----------------------------------------------------------------------
!> @brief
!> Write data to files
!-----------------------------------------------------------------------

MODULE WRITE_DATA

USE MPI

IMPLICIT NONE

CONTAINS

SUBROUTINE WRITE_ERROR(N, M, ER)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N, M     !< POLY ORDER
    INTEGER :: I, J
    
    DOUBLE PRECISION :: ER(0:N, 0:M) !< ERROR ON NODES
    
    OPEN(UNIT=1, FILE='error.dat')
    
    ! ERROR-------------------------------------------------------------
    DO J=0, M
        DO I=0, N
            WRITE(1, 10) I, J, ER(I, J)
        ENDDO
    ENDDO
    !-------------------------------------------------------------------
10 FORMAT(I2, 2X, I2, 2X, E20.10)
    CLOSE(UNIT=1)
    

END SUBROUTINE WRITE_ERROR

END MODULE WRITE_DATA
