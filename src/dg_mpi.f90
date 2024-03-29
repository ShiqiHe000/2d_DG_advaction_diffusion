!-----------------------------------------------------------------------
!> @brief
!> Activate OpenMPI
!-----------------------------------------------------------------------
MODULE SET_MPI

USE MPI
USE PARAM, ONLY: RANK, NUM_PROC, IERROR

IMPLICIT NONE

CONTAINS 

SUBROUTINE START_MPI

    ! MPI INITIALIZE 
    CALL MPI_INIT(IERROR)
    
    ! START UP MPI_ENVIRONEMENT
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERROR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUM_PROC, IERROR)

END SUBROUTINE START_MPI



END MODULE SET_MPI
