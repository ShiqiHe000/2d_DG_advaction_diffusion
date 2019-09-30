!-----------------------------------------------------------------------
!> @Generate the dual graph coordinates of the mesh file. n(i, j)
!! Coordinate expansion: 0 <= i <= EXP_X - 1, 0 <= j <= EXP_Y - 1
!-----------------------------------------------------------------------

MODULE GEN_DUAL_GRAPH

USE MPI
USE NODAL_2D_STORAGE
USE PARAM
USE GET_DUAL_COORD

IMPLICIT NONE 

CONTAINS

SUBROUTINE GEN_DUAL_GRAPH_2D

    IMPLICIT NONE 
    
    INTEGER :: K
    
    !-------------------------------------------------------------------
    NUM_OF_ELEMENT_X = 2**EXP_X
    NUM_OF_ELEMENT_Y = 2**EXP_Y
    NUM_OF_ELEMENT = NUM_OF_ELEMENT_X * NUM_OF_ELEMENT_Y
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    ALLOCATE(DUAL_COORD(2, NUM_OF_ELEMENT))
    
    DUAL_COORD = 0
    !-------------------------------------------------------------------
    
    ! ELEMENT SIZE START WITH UNIFORM-----------------------------------
    DELTA_X = (GX_R - GX_L) / DBLE(NUM_OF_ELEMENT_X)
    DELTA_Y = (GY_R - GY_L) / DBLE(NUM_OF_ELEMENT_Y)
    !-------------------------------------------------------------------
    
    DO K = 1, NUM_OF_ELEMENT
        CALL GET_DUAL_COORD_2D(ELEM_X_POSITION(1, K), &
                               ELEM_Y_POSITION(1, K), &
                               DELTA_X, DELTA_Y, DUAL_COORD(:, K))
    
    ENDDO
    

END SUBROUTINE GEN_DUAL_GRAPH_2D


END MODULE GEN_DUAL_GRAPH
