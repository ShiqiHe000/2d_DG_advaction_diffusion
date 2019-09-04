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
    USE NODAL_2D_STORAGE

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
    
    CALL GET_NODAL_2D_STORGAE_BASIS
    
    CALL GET_NODAL_2D_STORGAE_EXTENDS

    
    !-------------------------------------------------------------------
    
    ! END MPI-----------------------------------------------------------
    CALL MPI_FINALIZE(IERROR)
    !-------------------------------------------------------------------


END PROGRAM MAIN_LOOP
