!-----------------------------------------------------------------------
!> @brief 
!> The main loop of the program
!-----------------------------------------------------------------------

PROGRAM MAIN_LOOP
!-----------------------------------------------------------------------
!   POISSON PROBLEM
!-----------------------------------------------------------------------
    USE MPI
    USE SET_MPI
!    USE ADVECTION_DIFFUSION_DRIVER
!    USE VERIFICATION
    USE END_PROGRAM
    USE PREPARE_HILBERT_SCHEME
    use BASIS_STORAGE

    IMPLICIT NONE
    
    
    ! START MPI---------------------------------------------------------
    CALL START_MPI
    !-------------------------------------------------------------------
    
    ! PREPARE HILBERT CURVE---------------------------------------------
    CALL HILBERT_NUMBERING
    !-------------------------------------------------------------------
    
    ! CONOSTRUCT DG BASIS-----------------------------------------------
    CALL CONSTRUCT_BASIS_STORAGE
    !-------------------------------------------------------------------
    
    ! START GAMES-------------------------------------------------------
!    CALL DRIVER_FOR_DG_APPROXIMATION
    !-------------------------------------------------------------------
    
    ! VERIFICATION------------------------------------------------------
!    CALL GET_ERROR
    !-------------------------------------------------------------------
    
    ! FREE STORAGES-----------------------------------------------------
    CALL DEALLOCATE_ALL
    !-------------------------------------------------------------------
    
    
    ! END MPI-----------------------------------------------------------
    CALL MPI_FINALIZE(IERROR)
    !-------------------------------------------------------------------


END PROGRAM MAIN_LOOP
