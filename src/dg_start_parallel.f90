!-----------------------------------------------------------------------
!> @brief
!> Processor 1 finished reading mesh. 
!! Now to spread the load equally between processors (assume uniform
!! mesh and uniform polynomial order)
!-----------------------------------------------------------------------

MODULE DG_START_PARALLEL

USE MPI
USE FIRST_DISTRIBUTE
USE LOCAL_STORAGE
USE PARAM

IMPLICIT NONE 

CONTAINS

SUBROUTINE START_PARALLEL

!    INTEGER :: I

    ! DISTRIBUTE THE ELEMENTS EVENLY BETWEEN PROCESSORS-----------------
    CALL DISTRIBUTE_ELEM
    !-------------------------------------------------------------------
    
!    PRINT *, "RANK", RANK, "LOCAL_ELEM", LOCAL_ELEM_NUM
    
!    DO I=1, LOCAL_ELEM_NUM
!        PRINT *, X_LOCAL(:, I)
!    ENDDO


END SUBROUTINE START_PARALLEL

END MODULE DG_START_PARALLEL
