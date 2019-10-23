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

SUBROUTINE DEALLOCATE_NODAL_2D_STORAGE

    ! dg_nodal_2d_storage.f90-------------------------------------------
    DEALLOCATE(GL_POINT_X_T, GL_POINT_Y_T)
    DEALLOCATE(GL_W_X_T, GL_W_Y_T)
    
    DEALLOCATE(M_FIRST_DER_X_T, M_FIRST_DER_Y_T)
    
    DEALLOCATE(LAGRANGE_LEFT_T, LAGRANGE_RIGHT_T)
    DEALLOCATE(LAGRANGE_UP_T, LAGRANGE_DOWN_T)
    
    DEALLOCATE(SOLUTION)
    
    DEALLOCATE(PLEVEL_X, PLEVEL_Y)
    !-------------------------------------------------------------------


END SUBROUTINE DEALLOCATE_NODAL_2D_STORAGE

SUBROUTINE DEALLOCATE_ALL
    
    ! dg_nodal_2d_storage.f90-------------------------------------------
    DEALLOCATE(GL_POINT_X_T, GL_POINT_Y_T)
    DEALLOCATE(GL_W_X_T, GL_W_Y_T)
    
    DEALLOCATE(M_FIRST_DER_X_T, M_FIRST_DER_Y_T)
    
    DEALLOCATE(LAGRANGE_LEFT_T, LAGRANGE_RIGHT_T)
    DEALLOCATE(LAGRANGE_UP_T, LAGRANGE_DOWN_T)
    
    DEALLOCATE(SOLUTION)
    
    DEALLOCATE(PLEVEL_X, PLEVEL_Y)
    !-------------------------------------------------------------------
    
    ! dg_verification.f90-----------------------------------------------
!    DEALLOCATE(EXACT)
!    DEALLOCATE(ERROR)
!    DEALLOCATE(L2_NORM)
    !-------------------------------------------------------------------

END SUBROUTINE DEALLOCATE_ALL

END MODULE END_PROGRAM
