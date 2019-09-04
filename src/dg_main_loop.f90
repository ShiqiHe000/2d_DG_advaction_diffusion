!-----------------------------------------------------------------------
!> @brief 
!> The main loop of the program
!-----------------------------------------------------------------------

PROGRAM MAIN_LOOP
!-----------------------------------------------------------------------
!   POISSON PROBLEM
!-----------------------------------------------------------------------
    USE MPI
!    USE BASIS
    USE PARAM, ONLY: IERROR, MESHFILE
    USE SET_MPI
!    USE MESH
!    USE GRAPH_PARTITION
!    USE FIELDS
!    USE ADVECTION_DIFFUSION

    IMPLICIT NONE
    
!    INTEGER :: I
!    PRINT *, "CHECK"
    
    ! START MPI---------------------------------------------------------
    CALL START_MPI
    !-------------------------------------------------------------------
    
!    CALL READ_MESH(MESHFILE)
    
!    CALL BFS_NEW
    
!    CALL WRITE_TO_FILES
    !-------------------------------------------------------------------
    
!    PRINT *, "CHECK"
    
!    CALL DG_2D_WAVE
    
!    PRINT *, "CHECK"

    
    !-------------------------------------------------------------------
    
    ! END MPI-----------------------------------------------------------
    CALL MPI_FINALIZE(IERROR)
    !-------------------------------------------------------------------


END PROGRAM MAIN_LOOP
