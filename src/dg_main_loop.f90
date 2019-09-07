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
    USE ADVECTION_DIFFUSION_DRIVER
    USE VERIFICATION
    USE END_PROGRAM

    IMPLICIT NONE
    
!    INTEGER :: I
!    PRINT *, "CHECK"
    
    ! START MPI---------------------------------------------------------
    CALL START_MPI
    !-------------------------------------------------------------------
    
    ! START GAMES-------------------------------------------------------
    CALL DRIVER_FOR_DG_APPROXIMATION
    !-------------------------------------------------------------------
    
    ! VERIFICATION------------------------------------------------------
    CALL GET_ERROR
    !-------------------------------------------------------------------
    
    ! FREE STORAGES-----------------------------------------------------
    CALL DEALLOCATE_ALL
    !-------------------------------------------------------------------
    
    
    ! END MPI-----------------------------------------------------------
    CALL MPI_FINALIZE(IERROR)
    !-------------------------------------------------------------------


END PROGRAM MAIN_LOOP
