!-----------------------------------------------------------------------
!> @brief
!> Now let's end the games! Dellatocate all the variables.
!-----------------------------------------------------------------------

MODULE END_PROGRAM

USE MPI
USE NODAL_2D_STORAGE

IMPLICIT NONE

CONTAINS

SUBROUTINE DEALLOCATE_ALL

    DEALLOCATE(GL_POINT_X, GL_POINT_Y)
    DEALLOCATE(GL_W_X, GL_W_Y)
    
    DEALLOCATE(M_FIRST_DER_X, M_FIRST_DER_Y)
    
    DEALLOCATE(LAGRANGE_LEFT, LAGRANGE_RIGHT)
    DEALLOCATE(LAGRANGE_UP, LAGRANGE_DOWN)
    
    DEALLOCATE(SOLUTION)

END SUBROUTINE DEALLOCATE_ALL

END MODULE END_PROGRAM
