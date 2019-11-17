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
USE NODAL_2D_STORAGE

IMPLICIT NONE 

CONTAINS

SUBROUTINE START_PARALLEL

!    INTEGER :: I

    ! DISTRIBUTE THE ELEMENTS EVENLY BETWEEN PROCESSORS-----------------
    CALL DISTRIBUTE_ELEM
    !-------------------------------------------------------------------
    
    NUM_OF_ELEMENT_X = 2**EXP_X
    NUM_OF_ELEMENT_Y = 2**EXP_Y


END SUBROUTINE START_PARALLEL

END MODULE DG_START_PARALLEL
