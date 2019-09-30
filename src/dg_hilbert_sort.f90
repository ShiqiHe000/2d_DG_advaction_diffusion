!-----------------------------------------------------------------------
!> @brief
!> Re-numbering elements using Hilbert curve
!-----------------------------------------------------------------------

MODULE HILBERT_SORT

USE MPI
USE NODAL_2D_STORAGE
USE hilbert
USE PARAM, ONLY: EXP_X 

IMPLICIT NONE

CONTAINS

SUBROUTINE HIBERT_SORT_2D
    
    IMPLICIT NONE 
    
    INTEGER :: K
    INTEGER :: I, J     ! DUAL COORDINATES
    INTEGER :: D        ! HILBERT COORDINATE
    
    !-------------------------------------------------------------------
    ALLOCATE(X_HILBERT(4, NUM_OF_ELEMENT), Y_HILBERT(4, NUM_OF_ELEMENT))
    
    X_HILBERT = 0.0D0; Y_HILBERT = 0.0D0
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    DO K=1, NUM_OF_ELEMENT
    
        I = DUAL_COORD(1, K)
        J = DUAL_COORD(2, K)
        
        CALL xy2d ( EXP_X, I, J, D )
        
        X_HILBERT(:, D) = ELEM_X_POSITION(4, K)
        Y_HILBERT(:, D) = ELEM_Y_POSITION(4, K)
    
    ENDDO
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    DEALLOCATE(ELEM_X_POSITION, ELEM_Y_POSITION)
    DEALLOCATE(DUAL_COORD)
    !-------------------------------------------------------------------

END SUBROUTINE HIBERT_SORT_2D

END MODULE HILBERT_SORT
