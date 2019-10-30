!-----------------------------------------------------------------------
!> @brief
!> Now let's end the games! Dellatocate all the variables.
!-----------------------------------------------------------------------

MODULE END_PROGRAM

USE MPI
USE NODAL_2D_STORAGE
USE VERIFICATION
USE LOCAL_STORAGE

IMPLICIT NONE

CONTAINS

SUBROUTINE DEALLOCATE_NODAL_2D_STORAGE

    ! dg_nodal_2d_storage.f90-------------------------------------------
    DEALLOCATE(GL_POINT_X_T, GL_POINT_Y_T)
    DEALLOCATE(GL_W_X_T, GL_W_Y_T)
    
    DEALLOCATE(M_FIRST_DER_X_T, M_FIRST_DER_Y_T)
    
    DEALLOCATE(LAGRANGE_LEFT_T, LAGRANGE_RIGHT_T)
    DEALLOCATE(LAGRANGE_UP_T, LAGRANGE_DOWN_T)
    
    !-------------------------------------------------------------------


END SUBROUTINE DEALLOCATE_NODAL_2D_STORAGE




SUBROUTINE DEALLOCATE_LOCAL_STORAGE
    
    ! dg_local_storage.f90----------------------------------------------
    DEALLOCATE(SOLUTION)
    
    DEALLOCATE(PLEVEL_X, PLEVEL_Y)
    
    DEALLOCATE(A_LEVEL)
    
    DEALLOCATE(ELEM_RANGE)
    
    DEALLOCATE(X_LOCAL, Y_LOCAL)
    !------------------------------------------------------------------
    
END SUBROUTINE DEALLOCATE_LOCAL_STORAGE



SUBROUTINE DEALLOCATE_VERIFICATION

    DEALLOCATE(EXACT)
    DEALLOCATE(ERROR)
    DEALLOCATE(L2_NORM)

END SUBROUTINE DEALLOCATE_VERIFICATION



SUBROUTINE DEALLOCATE_ALL
    
    ! dg_nodal_2d_storage.f90-------------------------------------------
    CALL DEALLOCATE_NODAL_2D_STORAGE
    !-------------------------------------------------------------------
    
    ! dg_local_storage.f90----------------------------------------------
    CALL DEALLOCATE_LOCAL_STORAGE
    !-------------------------------------------------------------------

    
    ! dg_verification.f90-----------------------------------------------
    IF(VERIFICATION_SWITCH) CALL DEALLOCATE_VERIFICATION
    !-------------------------------------------------------------------

END SUBROUTINE DEALLOCATE_ALL



END MODULE END_PROGRAM
