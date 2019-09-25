!-----------------------------------------------------------------------
!> @brief
!> Construct mesh element by element
!-----------------------------------------------------------------------

MODULE GLOBAL_MESH

USE MPI
USE PARAM, ONLY: NUM_OF_ELEMENT_X, NUM_OF_ELEMENT_Y, NONE_L, NONE_R
USE DG_2D_CONSTRUCTOR
USE NODAL_2D_STORAGE, ONLY: ELEM_POINTER, NUM_OF_ELEMENT

IMPLICIT NONE

CONTAINS

SUBROUTINE GLOBAL_MESH_CONSTRUCT
!-----------------------------------------------------------------------
! Now all the element share the same basis information
! future work: element generates the basis bases on the local feature
!-----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: K
    
    !-------------------------------------------------------------------
    ALLOCATE(ELEM_POINTER(0:NUM_OF_ELEMENT, 2)) ! K+1 POINTS, LEFT AND RIGTH
    
    ELEM_POINTER = 0
    !-------------------------------------------------------------------
    
    CALL CONSTRUCT_BASIS
    
    ! LEFT BOUND IS ELEMENT K, RIGTH BOUND IS ELEMENT K+1
    DO K=1, NUM_OF_ELEMENT-1
        ELEM_POINTER(K, 1) = K
        ELEM_POINTER(K, 2) = K+1
    
    ENDDO
    
    ! FIRST AND LAST ELEMENT BOUNDARY
    ELEM_POINTER(0, 1) = NONE_L
    ELEM_POINTER(0, 2) = 1
    ELEM_POINTER(NUM_OF_ELEMENT, 1) = NUM_OF_ELEMENT
    ELEM_POINTER(NUM_OF_ELEMENT, 2) = NONE_R
    
    CALL INTERFACE_NODE_POSITION_CONSTRUCT


END SUBROUTINE GLOBAL_MESH_CONSTRUCT

END MODULE GLOBAL_MESH
