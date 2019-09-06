!-----------------------------------------------------------------------
!> @brief 
!> The main loop of the program
!-----------------------------------------------------------------------

PROGRAM MAIN_LOOP
!-----------------------------------------------------------------------
!   POISSON PROBLEM
!-----------------------------------------------------------------------
    USE MPI
    USE BASIS
    USE PARAM, ONLY: IERROR, MESHFILE
    USE SET_MPI
    USE DG_2D_CONSTRUCTOR

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
    
!    CALL CONSTRUCT_BASIS

    
    !-------------------------------------------------------------------
    
    ! END MPI-----------------------------------------------------------
    CALL MPI_FINALIZE(IERROR)
    !-------------------------------------------------------------------


END PROGRAM MAIN_LOOP
