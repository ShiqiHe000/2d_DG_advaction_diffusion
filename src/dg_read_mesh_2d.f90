!-----------------------------------------------------------------------
!> @brief
!> Read_mesh module read mesh information directly from GMSH .msh file.
!! NOTE: the .msh file has to be in version 2 ASII format.  
!-----------------------------------------------------------------------

MODULE READ_MESH

USE MPI

IMPLICIT NONE 

CONTAINS

SUBROUTINE READ_MESH_2D
!-----------------------------------------------------------------------
! READ .MSH FILE, RECORD THE XY COORDINATE OF EACH ELEMENT.
! SORTING THE NODE-ORDERING FORMAT. WARRENT EACH ELEMENT NODE IS ORDERING
! AS THE FORMAT:
! 4--------3
! |        |
! |        |
! |        |
! 1--------2
!
! GET ELEMENT IJ COORDINDATE.
! NOTE: THE ELEMENT SEQUENCE IN .MSH FILE DOES NOT MATTER.
!-----------------------------------------------------------------------
    
    USE PARAM, ONLY: MESHFILE, RANK
    
    IMPLICIT NONE 
    
    CHARACTER(LEN=50) :: CHARLINE   ! DUMMY LINE
    CHARACTER(LEN=50) :: LINE   ! STORE ELEMENT LINES
    
    INTEGER :: NUM_OF_PHY_NAME  ! NUMBER OF PHYSICAL NAME
    INTEGER :: TOTAL_NODE   ! TOTAL NUMBER OF NODES
    INTEGER :: NUM_OF_ELEMENT   ! NUMBER OF ELEMENTS
    INTEGER :: TOTAL_QUAD = 0     ! 4-NODE QUADRANGLE ELEMENT NUMBER
    
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: QUAD_NODE  ! QUADRANGLE NODES
    
    INTEGER :: TAG, ELEM_TYPE
    INTEGER :: I, A, B, C, D, E
    
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: NODE_XY ! NODE COORDINATES 
    
    IF(RANK == 0) THEN
    
        PRINT *, "------------------------------------------------------"
        PRINT *, "START READ MESH FILE"
        PRINT *, "------------------------------------------------------"
        
        OPEN(UNIT = 3, FILE = MESHFILE, STATUS = 'OLD')
        
        ! SKIP DUMMY LINES----------------------------------------------
        DO WHILE (.TRUE.)
          READ(1, *) CHARLINE
          CHARLINE = TRIM(CHARLINE)
          IF(CHARLINE == "$PhysicalNames") EXIT
        ENDDO
        !---------------------------------------------------------------
        
        READ(3, *) NUM_OF_PHY_NAME
        
        DO I = 1, NUM_OF_PHY_NAME
            READ(3, *) TAG
        ENDDO
        
        ! ERROR ALERT---------------------------------------------------
        IF (TAG /= 2) THEN
            PRINT *, "You forgot to set the physical surface."
            STOP
        ENDIF
        !---------------------------------------------------------------
        
        ! DUMMY LINES---------------------------------------------------
        READ(3, *)
        READ(3, *)
        !---------------------------------------------------------------
        
        READ(3, *) TOTAL_NODE
        
        ALLOCATE(NODE_XY(2, TOTAL_NODE))
        NODE_XY = 0.0D0
        
        ! READ NODE COORDINATES-----------------------------------------
        DO I=1, TOTAL_NODE
            READ(3, *) A, NODE_XY(1, I), NODE_XY(2, I)
        ENDDO
        !---------------------------------------------------------------
        
        !---------------------------------------------------------------
        READ(3, *)  ! READ DUMMY LINE, $EndNodes
        READ(3, *)  ! $Elements
        !---------------------------------------------------------------
        
        ! NUM OF ELEMENT------------------------------------------------
        READ(3, *) NUM_OF_ELEMENT
        !---------------------------------------------------------------
        
        ALLOCATE(QUAD_NODE(4, NUM_OF_ELEMENT))
        
        QUAD_NODE = 0
        
        ! READ ELEMENT NODES--------------------------------------------
        DO I=1, NUM_OF_ELEMENT
            READ(1, '(A)') LINE     ! STORE THE WHOLE LINE AS A STRING
                
            LINE = TRIM(LINE)
            
            READ(LINE, *) A, ELEM_TYPE  ! GET THE ELEMENT TYPE
            
            ! IF TYPE IS 4-NODE QUADRANGLE
            IF (ELEM_TYPE == 3) THEN
                TOTAL_QUAD = TOTAL_QUAD+1
                    
                READ(LINE, *) A, B, C, D, E, &
                            QUAD_NODE(1, TOTAL_QUAD), &
                            QUAD_NODE(2, TOTAL_QUAD), &
                            QUAD_NODE(3, TOTAL_QUAD), &
                            QUAD_NODE(4, TOTAL_QUAD)
            ENDIF 
            
            CALL SORT_NODE_ORDERING(NUM_OF_ELEMENT, TOTAL_QUAD, &
                                    QUAD_NODE, NODE_XY)
        
        ENDDO
        !---------------------------------------------------------------
        
        DEALLOCATE(QUAD_NODE)
        DEALLOCATE(NODE_XY)
        
    ENDIF
    

END SUBROUTINE READ_MESH_2D

SUBROUTINE SORT_NODE_ORDERING(TOTAL_NODE, NUM_OF_ELEMENTM, TOTAL_QUAD, &
                                QUAD_NODE, NODE_XY)
!-----------------------------------------------------------------------
! SORTING THE NODE-ORDERING FORMAT. WARRENT EACH ELEMENT NODE IS ORDERING
! AS THE FORMAT:
! 4--------3
! |        |
! |        |
! |        |
! 1--------2
!-----------------------------------------------------------------------
    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: NUM_OF_ELEMENT
    INTEGER, INTENT(IN) :: TOTAL_QUAD
    INTEGER, INTENT(IN) :: TOTAL_NODE
    
    INTEGER :: QUAD_NODE(4, NUM_OF_ELEMENT)
    
    INTEGER :: K, I
    
    DOUBLE PRECISION :: NODE_XY(2, TOTAL_NODE)
    
    DO K=1, TOTAL_QUAD
        
    
    ENDDO
    

END SUBROUTINE SORT_NODE_ORDERING


END MODULE READ_MESH
