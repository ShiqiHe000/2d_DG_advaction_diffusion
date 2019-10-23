!-----------------------------------------------------------------------
!> @brief
!> Processor 1 finished reading mesh. 
!! Now to spread the load equally between processors (assume uniform
!! mesh and uniform polynomial order)
!-----------------------------------------------------------------------

MODULE DG_START_PARALLEL

USE MPI
USE FIRST_DISTRIBUTE

IMPLICIT NONE 

CONTAINS

SUBROUTINE START_PARALLEL

    ! DISTRIBUTE THE ELEMENTS EVENLY BETWEEN PROCESSORS-----------------
    CALL DISTRIBUTE_ELEM
    !-------------------------------------------------------------------
    
    

END SUBROUTINE START_PARALLEL

END MODULE DG_START_PARALLEL
