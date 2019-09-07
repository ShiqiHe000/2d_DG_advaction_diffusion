!-----------------------------------------------------------------------
!> @brief
!> Now let's end the games! Dellatocate all the variables.
!-----------------------------------------------------------------------

MODULE END_PROGRAM

USE MPI
USE NODAL_2D_STORAGE
USE VERIFICATION

IMPLICIT NONE

CONTAINS

SUBROUTINE DEALLOCATE_ALL
    
    ! dg_nodal_2d_storage.f90-------------------------------------------
    DEALLOCATE(GL_POINT_X, GL_POINT_Y)
    DEALLOCATE(GL_W_X, GL_W_Y)
    
    DEALLOCATE(M_FIRST_DER_X, M_FIRST_DER_Y)
    
    DEALLOCATE(LAGRANGE_LEFT, LAGRANGE_RIGHT)
    DEALLOCATE(LAGRANGE_UP, LAGRANGE_DOWN)
    
    DEALLOCATE(SOLUTION)
    !-------------------------------------------------------------------
    
    ! dg_verification.f90-----------------------------------------------
    DEALLOCATE(EXACT)
    DEALLOCATE(ERROR)
    DEALLOCATE(L2_NORM)
    !-------------------------------------------------------------------

END SUBROUTINE DEALLOCATE_ALL

END MODULE END_PROGRAM
