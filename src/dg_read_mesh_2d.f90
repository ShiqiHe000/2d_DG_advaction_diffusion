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
    
    USE PARAM, ONLY: MESHFILE, RANK, FILEPLACE
    
    IMPLICIT NONE 
    
    CHARACTER(LEN=50) :: CHARLINE   ! DUMMY LINE
    CHARACTER(LEN=50) :: LINE   ! STORE ELEMENT LINES
    
    INTEGER :: NUM_OF_PHY_NAME  ! NUMBER OF PHYSICAL NAME
    INTEGER :: TOTAL_NODE   ! TOTAL NUMBER OF NODES
    INTEGER :: NUM_OF_ELEMENT   ! NUMBER OF ELEMENTS
    INTEGER :: TOTAL_QUAD      ! 4-NODE QUADRANGLE ELEMENT NUMBER
    
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: QUAD_NODE  ! QUADRANGLE NODES
    
    INTEGER :: TAG, ELEM_TYPE
    INTEGER :: I, A, B, C, D, E
    
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: NODE_XY ! NODE COORDINATES 
    
    IF(RANK == 0) THEN
    
        PRINT *, "------------------------------------------------------"
        PRINT *, "START READ MESH FILE"
        PRINT *, "------------------------------------------------------"
        
        OPEN(UNIT = 3, FILE = FILEPLACE//MESHFILE, STATUS = 'OLD')
        
        ! SKIP DUMMY LINES----------------------------------------------
        DO WHILE (.TRUE.)
          READ(3, *) CHARLINE
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
        TOTAL_QUAD = 0
        
        ! READ ELEMENT NODES--------------------------------------------
        DO I=1, NUM_OF_ELEMENT
            READ(3, '(A)') LINE     ! STORE THE WHOLE LINE AS A STRING
                
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
            
        ENDDO
        !---------------------------------------------------------------
        
        CALL SORT_NODE_ORDERING(TOTAL_NODE, TOTAL_QUAD, &
                                    QUAD_NODE(:, 1:TOTAL_QUAD), NODE_XY)

        DEALLOCATE(QUAD_NODE)
        DEALLOCATE(NODE_XY)

        
    ENDIF
    

END SUBROUTINE READ_MESH_2D

SUBROUTINE SORT_NODE_ORDERING(TOTAL_NODE, TOTAL_QUAD, &
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
    USE NODAL_2D_STORAGE, ONLY: ELEM_X_POSITION, ELEM_Y_POSITION
!    USE WRITE_DATA
    
    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: TOTAL_QUAD
    INTEGER, INTENT(IN) :: TOTAL_NODE
    
    INTEGER :: X_SCORES(2), Y_SCORES(2) ! SCORES FOR THE RELATIVE POSITION
    
    INTEGER :: NODE1, NODE2, NODE3, NODE4   ! AIMED NODE SEQUENCES
    
    INTEGER :: QUAD_NODE(4, TOTAL_QUAD) ! NOTE: ONLY INPUT THE TOTAL_QUADS ELEMENT, TOTAL_QUAD < NUM_OF_ELEMENT
    
    INTEGER :: K, I, SCORE 
    
    DOUBLE PRECISION :: NODE_XY(2, TOTAL_NODE)
    
    DOUBLE PRECISION :: X_MAX, X_MIN    ! X COORDINATES
    DOUBLE PRECISION :: Y_MAX, Y_MIN    ! Y COORDINATES
    
    !-------------------------------------------------------------------
    X_SCORES = (/2, 1/)
    Y_SCORES = (/20, 10/)
    
    NODE1 = X_SCORES(2) + Y_SCORES(2)
    NODE2 = X_SCORES(2) + Y_SCORES(1)
    NODE3 = X_SCORES(1) + Y_SCORES(1)
    NODE4 = X_SCORES(1) + Y_SCORES(2)
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    ALLOCATE(ELEM_X_POSITION(4, TOTAL_QUAD))
    ALLOCATE(ELEM_Y_POSITION(4, TOTAL_QUAD))
    
    ELEM_X_POSITION = 0.0D0; ELEM_Y_POSITION = 0.0D0
    !-------------------------------------------------------------------
    
    DO K=1, TOTAL_QUAD
        CALL GET_STANDARD(NODE_XY(:, QUAD_NODE(1, K)), &
                          NODE_XY(:, QUAD_NODE(3, K)), &
                          X_MAX, X_MIN, Y_MAX, Y_MIN)
                          
        DO I=1, 4
            CALL GET_SCORES(X_MAX, Y_MAX, &
                            NODE_XY(:, QUAD_NODE(I, K)), &
                            SCORE, X_SCORES, Y_SCORES)
        
            IF(SCORE == NODE1) THEN
                ELEM_X_POSITION(1, K) = NODE_XY(1, QUAD_NODE(I, K))
                ELEM_Y_POSITION(1, K) = NODE_XY(2, QUAD_NODE(I, K))
            ELSEIF(SCORE == NODE2) THEN
                ELEM_X_POSITION(2, K) = NODE_XY(1, QUAD_NODE(I, K))
                ELEM_Y_POSITION(2, K) = NODE_XY(2, QUAD_NODE(I, K))
                
            ELSEIF(SCORE == NODE3) THEN
                ELEM_X_POSITION(3, K) = NODE_XY(1, QUAD_NODE(I, K))
                ELEM_Y_POSITION(3, K) = NODE_XY(2, QUAD_NODE(I, K))
                
            ELSEIF(SCORE == NODE4) THEN
                ELEM_X_POSITION(4, K) = NODE_XY(1, QUAD_NODE(I, K))
                ELEM_Y_POSITION(4, K) = NODE_XY(2, QUAD_NODE(I, K))
                
            ELSE
                PRINT *, "In dg_read_mesh_2d.f90, &
                            & the 'score' in SORT_NODE_ORDERING &
                            & does not match with the standards."
                STOP
            ENDIF
        
        ENDDO
        
        
        
    ENDDO
    
!    CALL WRITE_VISUAL(TOTAL_QUAD, ELEM_X_POSITION, ELEM_Y_POSITION)
    

END SUBROUTINE SORT_NODE_ORDERING

SUBROUTINE GET_STANDARD(NODE_XY1, NODE_XY3, X_MAX, X_MIN, Y_MAX, Y_MIN)
    
    IMPLICIT NONE 
    
    DOUBLE PRECISION :: X_MAX, X_MIN    ! X COORDINATES
    DOUBLE PRECISION :: Y_MAX, Y_MIN    ! Y COORDINATES
    
    DOUBLE PRECISION :: NODE_XY1(2)     ! FIRST NODE COORDINATE
    DOUBLE PRECISION :: NODE_XY3(2)     ! THIRD NODE COORDINATE
    
    DOUBLE PRECISION :: X1, Y1, X3, Y3  ! DIAGONAL NODES COORDINATES
    
    X1 = NODE_XY1(1); X3 = NODE_XY3(1)
    Y1 = NODE_XY1(2); Y3 = NODE_XY3(2)
    
    !-------------------------------------------------------------------
    IF(X1 > X3) THEN
        X_MAX = X1
        X_MIN = X3
    ELSE
        X_MAX = X3
        X_MIN = X1
    ENDIF
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    IF(Y1 > Y3) THEN
        Y_MAX = Y1
        Y_MIN = Y3
    ELSE
        Y_MAX = Y3
        Y_MIN = Y1
    ENDIF
    !-------------------------------------------------------------------
    
    
END SUBROUTINE GET_STANDARD

SUBROUTINE GET_SCORES(X_MAX, Y_MAX, NODE_XY1, SCORE, X_SCORES, Y_SCORES)

    IMPLICIT NONE 

    INTEGER :: SCORE 
    
    INTEGER, INTENT(IN) :: X_SCORES(2), Y_SCORES(2) ! SCORES FOR THE RELATIVE POSITION
    
    DOUBLE PRECISION :: X_MAX, Y_MAX
    
    DOUBLE PRECISION :: NODE_XY1(2)    
    
    LOGICAL :: FLAG1, FLAG2 !< LOGICAL OPERATOR
    
    SCORE = 0
    
    CALL ALMOSTEQUAL_2(FLAG1, NODE_XY1(1), X_MAX)
    CALL ALMOSTEQUAL_2(FLAG2, NODE_XY1(2), Y_MAX)

    IF(FLAG1) THEN
        SCORE = SCORE + X_SCORES(1)
    ELSE
        SCORE = SCORE + X_SCORES(2)
    ENDIF
    
    IF(FLAG2) THEN
        SCORE = SCORE + Y_SCORES(1)
    ELSE
        SCORE = SCORE + Y_SCORES(2)
    ENDIF
    


END SUBROUTINE GET_SCORES

SUBROUTINE ALMOSTEQUAL_2(FLAG, A, B)
!-----------------------------------------------------------------------
! DETERMINE WHETHER TWO FLOATING POINT NUMBER A AND B ARE EQUAL OR NOT
! RETURN LOGICAL FLAG
!-----------------------------------------------------------------------

    IMPLICIT NONE 
    
    DOUBLE PRECISION, INTENT(IN) :: A, B    ! TWO TARGETS TO COMPARE
    DOUBLE PRECISION, PARAMETER :: THRESHOLD = 1.0E-11
    
    LOGICAL :: FLAG 
    
    
    IF ( DABS(A - B) <= THRESHOLD ) THEN
        FLAG = .TRUE.
    ELSE 
        FLAG = .FALSE.
    ENDIF


END SUBROUTINE ALMOSTEQUAL_2

END MODULE READ_MESH
