!-----------------------------------------------------------------------
!> @brief
!> Write data to files
!-----------------------------------------------------------------------

MODULE WRITE_DATA

USE MPI
USE NODAL_2D_STORAGE
USE POLY_LEVEL_AND_ORDER
USE PARAM
USE POLY_LEVEL_AND_ORDER

IMPLICIT NONE

INTEGER :: FRAME = 1

CONTAINS

SUBROUTINE WRITE_ERROR(ER)

    IMPLICIT NONE
    
    INTEGER :: I, J, K
    INTEGER :: PORDERX, PORDERY
    
    DOUBLE PRECISION, INTENT(IN) :: ER(0:NMAX, 0:MMAX, NUM_OF_EQUATION, 0:NUM_OF_ELEMENT-1) !< ERROR ON NODES
    
    OPEN(UNIT=1, FILE='error.dat')
    
    ! ERROR-------------------------------------------------------------
    DO K = 0, NUM_OF_ELEMENT-1
        CALL POLY_LEVEL_TO_ORDER(N, PLEVEL_X(K), PORDERX)
        CALL POLY_LEVEL_TO_ORDER(M, PLEVEL_Y(K), PORDERY)
    
        DO J=0, PORDERY
            DO I=0, PORDERX
                WRITE(1, 10) I, J, ER(I, J, 1, K), ER(I, J, 2, K), ER(I, J, 3, K)
            ENDDO
        ENDDO
        
    ENDDO
    !-------------------------------------------------------------------
10 FORMAT(I2, 2X, I2, 2X, E20.10, 2X, E20.10, 2X, E20.10)
    CLOSE(UNIT=1)
    

END SUBROUTINE WRITE_ERROR

SUBROUTINE WRITE_RESULTS(N, M, EX, RESULTS)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N, M     !< POLY ORDER
    INTEGER :: I, J
    
    DOUBLE PRECISION :: EX(0:N, 0:M)    !< EXACT SOLUTION
    DOUBLE PRECISION :: RESULTS(0:N, 0:M)   !< COMPUTED RESULTS
    
    OPEN(UNIT=2, FILE='solutions.dat')
    
    DO J=0, M
        DO I=0, N
            WRITE(2, 20) I, J, RESULTS(I, J), EX(I, J)
        ENDDO
    ENDDO
20 FORMAT(I2, 2X, I2, 2X, E20.10, 2X, E20.10)
    CLOSE(UNIT=2)
END SUBROUTINE WRITE_RESULTS

!-----------------------------------------------------------------------
! Element number start with 0
!-----------------------------------------------------------------------
SUBROUTINE WRITE_VISUAL2(NEL_TOTAL, X_GLOBAL, Y_GLOBAL)
    
    USE PARAM, ONLY: RANK

    IMPLICIT NONE 
    
    CHARACTER(LEN=16) :: FILENAME
    
    INTEGER :: ELEM, IEL
    INTEGER, INTENT(IN) :: NEL_TOTAL    !< TOTAL NUMBER OF ELEMENT
    
    DOUBLE PRECISION :: X_GLOBAL(4, 0:NEL_TOTAL-1)  !< ELEMENT X COORDINATES
    DOUBLE PRECISION :: Y_GLOBAL(4, 0:NEL_TOTAL-1)  !< ELEMENT Y COORDINATES
    
    IF (RANK == 0) THEN
    
        WRITE(FILENAME,FMT='(''aoutput'',I5.5,''.dat'')') FRAME
        
        OPEN(4,FILE=FILENAME)
        
        WRITE(4, FMT='(''TITLE = "MESH AND SOLUTIONS"'')')
        WRITE(4, FMT='(''VARIABLES = "X", "Y"'')')
    
        ELEM=1
        
        DO IEL=0, NEL_TOTAL-1

            WRITE(4, 30) "ZONE T= ", '"', "IEL", ELEM, '"', &
                        "I=2, J=2","DATAPACKING = POINT" 
30 FORMAT(A8, A1, A3, I6, A1, 2X, A8, 2X, A19)

            ELEM = ELEM+1
            
            WRITE(4, 40) X_GLOBAL(1, IEL), Y_GLOBAL(1, IEL)
            WRITE(4, 40) X_GLOBAL(2, IEL), Y_GLOBAL(2, IEL)
            WRITE(4, 40) X_GLOBAL(4, IEL), Y_GLOBAL(4, IEL)
            WRITE(4, 40) X_GLOBAL(3, IEL), Y_GLOBAL(3, IEL)
            
            
        ENDDO
40 FORMAT(F10.5, 2X, F10.5, 2X, I2)

    CLOSE(UNIT=4)
    ENDIF


END SUBROUTINE WRITE_VISUAL2

!-----------------------------------------------------------------------
! Element number start with 1
!-----------------------------------------------------------------------
SUBROUTINE WRITE_VISUAL(NEL_TOTAL, X_GLOBAL, Y_GLOBAL)
    
    USE PARAM, ONLY: RANK

    IMPLICIT NONE 
    
    CHARACTER(LEN=16) :: FILENAME
    
    INTEGER :: ELEM, IEL
    INTEGER, INTENT(IN) :: NEL_TOTAL    !< TOTAL NUMBER OF ELEMENT
    
    DOUBLE PRECISION :: X_GLOBAL(4, NEL_TOTAL)  !< ELEMENT X COORDINATES
    DOUBLE PRECISION :: Y_GLOBAL(4, NEL_TOTAL)  !< ELEMENT Y COORDINATES
    
    IF (RANK == 0) THEN
    
        WRITE(FILENAME,FMT='(''aoutput'',I5.5,''.dat'')') FRAME
        
        OPEN(4,FILE=FILENAME)
        
        WRITE(4, FMT='(''TITLE = "MESH AND SOLUTIONS"'')')
        WRITE(4, FMT='(''VARIABLES = "X", "Y"'')')
    
        ELEM=1
        
        DO IEL=1, NEL_TOTAL

            WRITE(4, 30) "ZONE T= ", '"', "IEL", ELEM, '"', &
                        "I=2, J=2","DATAPACKING = POINT" 
30 FORMAT(A8, A1, A3, I6, A1, 2X, A8, 2X, A19)

            ELEM = ELEM+1
            
            WRITE(4, 40) X_GLOBAL(1, IEL), Y_GLOBAL(1, IEL)
            WRITE(4, 40) X_GLOBAL(2, IEL), Y_GLOBAL(2, IEL)
            WRITE(4, 40) X_GLOBAL(4, IEL), Y_GLOBAL(4, IEL)
            WRITE(4, 40) X_GLOBAL(3, IEL), Y_GLOBAL(3, IEL)
            
            
        ENDDO
40 FORMAT(F10.5, 2X, F10.5, 2X, I2)

    CLOSE(UNIT=4)
    ENDIF


END SUBROUTINE WRITE_VISUAL

END MODULE WRITE_DATA
