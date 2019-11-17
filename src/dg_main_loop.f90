!-----------------------------------------------------------------------
!> @brief 
!> The main loop of the program
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!   2D linear advaction problem
!-----------------------------------------------------------------------
PROGRAM MAIN_LOOP
    USE MPI
    USE SET_MPI
    USE ADVECTION_DIFFUSION_DRIVER
    USE VERIFICATION
    USE END_PROGRAM
    USE PREPARE_HILBERT_SCHEME
    USE PARAM
    USE DG_START_PARALLEL

    IMPLICIT NONE
    
    
    ! START MPI---------------------------------------------------------
    CALL START_MPI
    !-------------------------------------------------------------------
    
    ! PREPARE HILBERT CURVE---------------------------------------------
    CALL HILBERT_NUMBERING
    !-------------------------------------------------------------------
    
    ! START PARALLEL PROCESS -------------------------------------------
    CALL START_PARALLEL
    !-------------------------------------------------------------------
    
    ! START GAMES-------------------------------------------------------
    CALL DRIVER_FOR_DG_APPROXIMATION
    !-------------------------------------------------------------------
    
    ! VERIFICATION------------------------------------------------------
    IF (VERIFICATION_SWITCH) THEN
        CALL GET_ERROR
    ENDIF
    !-------------------------------------------------------------------
    
    ! FREE STORAGES-----------------------------------------------------
    CALL DEALLOCATE_ALL
    !-------------------------------------------------------------------
    
    
    ! END MPI-----------------------------------------------------------
    CALL MPI_FINALIZE(IERROR)
    !-------------------------------------------------------------------


END PROGRAM MAIN_LOOP
