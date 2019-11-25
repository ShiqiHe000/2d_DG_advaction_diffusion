!-----------------------------------------------------------------------
!> @brief
!> Output solution to files.
!! Serial version.
!-----------------------------------------------------------------------

MODULE OUTPUT

USE MPI
USE PARAM, ONLY: NUM_OF_EQUATION, N, M, RANK, MMAX, NMAX, RANK
USE POLY_LEVEL_AND_ORDER
USE INTERFACES_CONSTRUCT
USE NODAL_2D_STORAGE, ONLY: LAGRANGE_LEFT_T, LAGRANGE_RIGHT_T, &
                                LAGRANGE_DOWN_T, LAGRANGE_UP_T, &
                                GL_POINT_X_T, GL_POINT_Y_T
USE LOCAL_STORAGE

IMPLICIT NONE 

INTEGER :: FILE_NUM = 1 ! FILE NUMBER

CONTAINS

!-----------------------------------------------------------------------
!> OUTPUT MESHES AND SOLUTION TO TECPLOT FILES
!-----------------------------------------------------------------------
SUBROUTINE WRITE_MESH(NEL_TOTAL, X_GLOBAL, Y_GLOBAL, PLEVELX, PLEVELY, &
                        SOLUTION_ALL, T)
                        
    USE PARAM, ONLY: OUTPUT_PLACE
    
    IMPLICIT NONE 
    
    CHARACTER(LEN=16) :: FILENAME
    
    INTEGER :: ELEM, IEL
    INTEGER :: PORDERX, PORDERY
    INTEGER, INTENT(IN) :: NEL_TOTAL    !< TOTAL NUMBER OF ELEMENT
    
    INTEGER :: PLEVELX(0:NEL_TOTAL-1)  !< POLYNOMIAL LEVEL IN X(CAN BE TRANSFORM TO POLY ORDER)
    INTEGER :: PLEVELY(0:NEL_TOTAL-1)  !< POLYNOMIAL LEVEL IN Y(CAN BE TRANSFORM TO POLY ORDER)
    
    
    DOUBLE PRECISION :: X_GLOBAL(2, 0:NEL_TOTAL-1)  !< ELEMENT X COORDINATES
    DOUBLE PRECISION :: Y_GLOBAL(2, 0:NEL_TOTAL-1)  !< ELEMENT Y COORDINATES
    
    DOUBLE PRECISION :: SOLUTION_ALL(0:NMAX, 0:MMAX, NUM_OF_EQUATION, 0:NEL_TOTAL-1)    !< NUMERICAL APPROXIMATION
    
    DOUBLE PRECISION :: T !< CURRENT TIME
    
    DOUBLE PRECISION :: SOLU_INT_L(0:MMAX, NUM_OF_EQUATION)
    DOUBLE PRECISION :: SOLU_INT_R(0:MMAX, NUM_OF_EQUATION)
    
    DOUBLE PRECISION :: SOLU_INT_D(0:NMAX, NUM_OF_EQUATION)
    DOUBLE PRECISION :: SOLU_INT_U(0:NMAX, NUM_OF_EQUATION)
    
    
    SOLU_INT_L = 0.0D0; SOLU_INT_R = 0.0D0
    SOLU_INT_D = 0.0D0; SOLU_INT_U = 0.0D0
    
    ! OPEN FILE
    WRITE(FILENAME,FMT='(''aoutput'',I5.5,''.dat'')') FILE_NUM
    
    ! PROC 1 WRITES HEADER----------------------------------------------
    IF (RANK == 0) THEN
    
        OPEN(5,FILE = OUTPUT_PLACE//FILENAME)
        
        WRITE(5, FMT='(''TITLE = "MESH AND SOLUTIONS"'')')
        WRITE(5, FMT='(''VARIABLES = "X", "Y", "PRESSURE", "U", "V", &
                        & "N", "M", "PROC"'')')
                        
        ELEM = 1
    ELSE
    
        OPEN(5,FILE = OUTPUT_PLACE//FILENAME, ACCESS = 'APPEND', STATUS = 'OLD')
        ELEM = ELEM_RANGE(RANK) + 1
    ENDIF
    !-------------------------------------------------------------------
    
    
    DO IEL=0, NEL_TOTAL-1
    
        CALL POLY_LEVEL_TO_ORDER(N, PLEVELX(IEL), PORDERX)
        CALL POLY_LEVEL_TO_ORDER(M, PLEVELY(IEL), PORDERY)

        WRITE(5, 50) "ZONE T= ", '"', "IEL", ELEM, '",', &
                    "I=", PORDERX+3, ',', "J=", PORDERY+3, ',', &
                    "SOLUTIONTIME=", T, ',', &
                    "DATAPACKING = POINT" 
50 FORMAT(A8, A1, A3, I6, A2, 2X, A2, I2, A1, 2X, A2, I2, A1, 2X, A13, F10.5, A1, 2X, A19)
            
        ELEM = ELEM+1
        
        ! INTERPOLATE TO INTERFACES---------------------------------
        ! X DIRECTION
        CALL CONSTRUCT_INTERFACES_X(PORDERX, PORDERY, &
                                        NUM_OF_EQUATION, &
                                        SOLUTION_ALL(0:PORDERX, 0:PORDERY, :, IEL), &
                                        LAGRANGE_LEFT_T(0:PORDERX, PLEVELX(IEL)),&
                                        LAGRANGE_RIGHT_T(0:PORDERX, PLEVELX(IEL)), &
                                        SOLU_INT_L(0:PORDERY, :), &
                                        SOLU_INT_R(0:PORDERY, :) )
        
        ! Y DIRECTION
        CALL CONSTRUCT_INTERFACES_Y(PORDERX, PORDERY, &
                                        NUM_OF_EQUATION, &
                                        SOLUTION_ALL(0:PORDERX, 0:PORDERY, :, IEL), &
                                        LAGRANGE_DOWN_T(0:PORDERY, PLEVELY(IEL)),&
                                        LAGRANGE_UP_T(0:PORDERY, PLEVELY(IEL)), &
                                        SOLU_INT_D(0:PORDERX, :), &
                                        SOLU_INT_U(0:PORDERX, :) )
        !-----------------------------------------------------------
        
        !-----------------------------------------------------------
        ! X INTERFACE (LEFT)
        CALL X_INTERFACE(PORDERX, PORDERY, PLEVELX(IEL), &
                            SOLU_INT_D(0:PORDERX, :), &
                            X_GLOBAL(1, IEL), X_GLOBAL(2, IEL), &
                            Y_GLOBAL(1, IEL))
        
        ! INTERIOR NODES (INCLUDING NODES ON Y INTERFACES)
        CALL INTERIOR_MESH(PORDERX, PORDERY, PLEVELX(IEL), &
                            PLEVELY(IEL), &
                            SOLUTION_ALL(0:PORDERX, 0:PORDERY, :, IEL), &
                            SOLU_INT_L(0:PORDERY, :), SOLU_INT_R(0:PORDERY, :), &
                            X_GLOBAL(1, IEL), X_GLOBAL(2, IEL), &
                            Y_GLOBAL(1, IEL), Y_GLOBAL(2, IEL))
        
        ! X INTERFACE (RIGHT)
        CALL X_INTERFACE(PORDERX, PORDERY, PLEVELX(IEL), &
                            SOLU_INT_U(0:PORDERX, :), &
                            X_GLOBAL(1, IEL), X_GLOBAL(2, IEL), &
                            Y_GLOBAL(2, IEL))
        !-----------------------------------------------------------
    
    ENDDO
    
    FILE_NUM = FILE_NUM + 1
    CLOSE(UNIT=5)


END SUBROUTINE WRITE_MESH

!-----------------------------------------------------------------------
!> INTERPOLATE THE SOLUTIONS ON INTERFACES TO THE CORNERS.
!! NOTE: ONLY SUPPORT INTERFACES ON X DIRECTION.
!-----------------------------------------------------------------------
SUBROUTINE INTERPOLATE_TO_CORNER(N_NOW, PLEVEL, Q, CORNER)

    USE BASIS
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N_NOW    !< INTERFACE POLY ORDER
    INTEGER, INTENT(IN) :: PLEVEL    !< POLY LEVEL
    
    INTEGER :: S
    
    DOUBLE PRECISION :: Q(0:N_NOW, NUM_OF_EQUATION) !< SOLUTIONS ON THE INTERFACE COLLOCATION NODES 
    
    DOUBLE PRECISION :: CORNER(2, NUM_OF_EQUATION)   !< INTERPOTATES TO TWO CORNERS
    
    CORNER = 0.0D0
    
    DO S = 1, NUM_OF_EQUATION
        CALL INTERPOLATE_TO_BOUNDARY(N_NOW, Q(:, S), &
                                    LAGRANGE_LEFT_T(0:N_NOW, PLEVEL), &
                                    CORNER(1, S))
        CALL INTERPOLATE_TO_BOUNDARY(N_NOW, Q(:, S), &
                                    LAGRANGE_RIGHT_T(0:N_NOW, PLEVEL), &
                                    CORNER(2, S))
    ENDDO


END SUBROUTINE INTERPOLATE_TO_CORNER

!-----------------------------------------------------------------------
!> CONSTRUCT MESH ON X INTERFACES
!-----------------------------------------------------------------------
SUBROUTINE X_INTERFACE(N_NOW, M_NOW, PLEVEL, Q, X1, X2, Y_STALL)

     USE AFFINE_MAP
    
    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: N_NOW, M_NOW    !< INTERFACE POLY ORDER
    INTEGER, INTENT(IN) :: PLEVEL    !< POLY LEVEL    
    INTEGER :: I, S
    
    DOUBLE PRECISION :: Q(0:N_NOW, NUM_OF_EQUATION) !< SOLUTIONS ON THE INTERFACE COLLOCATION NODES 
    
    DOUBLE PRECISION :: CORNER(2, NUM_OF_EQUATION)  !< SOLUTION ON THE CORNERS
    
    DOUBLE PRECISION :: X1, X2  !< X COORDINATES ON TWO CORNERS(PHYSICAL)
    DOUBLE PRECISION :: Y_STALL !< Y COORDINATES (PHYSICAL, DOES NOT CHANGE)
    
    DOUBLE PRECISION :: DEL_X   ! ELEMENT SIZE IN X DIRECTION
    DOUBLE PRECISION :: X
    
    DEL_X = X2 - X1
    
    CALL INTERPOLATE_TO_CORNER(N_NOW, PLEVEL, Q, CORNER)
    
    ! CORNER 1
    WRITE(5, FMT = '(F10.5, 2X, F10.5, 2X)', ADVANCE = 'NO') X1, Y_STALL 
    DO S = 1, NUM_OF_EQUATION
        WRITE(5, FMT = '(F10.5, 2X)', ADVANCE = 'NO') CORNER(1, S)
    ENDDO
    WRITE(5, FMT = '(3I6)') N_NOW, M_NOW, RANK
    
    ! COLLOCATION POINTS
    DO I = 0, N_NOW
        CALL AFFINE_MAPPING(GL_POINT_X_T(I, PLEVEL), &
                            X, X1, DEL_X)
                            
        WRITE(5, FMT = '(F10.5, 2X, F10.5, 2X)', ADVANCE = 'NO') X, Y_STALL
        
        DO S = 1, NUM_OF_EQUATION
            WRITE(5, FMT = '(F10.5, 2X)', ADVANCE = 'NO') Q(I, S)
        ENDDO
        
        WRITE(5, FMT = '(3I6)') N_NOW, M_NOW, RANK
    ENDDO
    
    ! CORNER 2
    WRITE(5, FMT = '(F10.5, 2X, F10.5, 2X)', ADVANCE = 'NO') X2, Y_STALL 
    DO S = 1, NUM_OF_EQUATION
        WRITE(5, FMT = '(F10.5, 2X)', ADVANCE = 'NO') CORNER(2, S)
    ENDDO
    WRITE(5, FMT = '(3I6)') N_NOW, M_NOW, RANK
    
END SUBROUTINE X_INTERFACE

SUBROUTINE INTERIOR_MESH(N_NOW, M_NOW, PLX, PLY, Q, Q_Y1, Q_Y2, &
                            X1, X2, Y1, Y2)
    USE AFFINE_MAP
    
    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: N_NOW, M_NOW !< POLY ORDER
    INTEGER, INTENT(IN) :: PLX    !< POLY LEVEL IN X DIRECTION
    INTEGER, INTENT(IN) :: PLY    !< POLY LEVEL IN Y DIRECTION
    
    INTEGER :: I, J, S
    
    DOUBLE PRECISION :: Q(0:N_NOW, 0:M_NOW, NUM_OF_EQUATION)    !< SOLUTION ON INTERIOR NODES(INCLUDING NODES ON Y INTERFACES)
    DOUBLE PRECISION :: Q_Y1(0:M_NOW, NUM_OF_EQUATION)  !< SOLUTION ON Y INTERFACE (BOTTOM)
    DOUBLE PRECISION :: Q_Y2(0:M_NOW, NUM_OF_EQUATION)  !< SOLUTION ON Y INTERFACE (TOP)
    
    DOUBLE PRECISION, INTENT(IN) :: X1, Y1, X2, Y2  !< NODE1 XY COORDS OF CURRENT ELEMENT
    DOUBLE PRECISION :: X, Y
    
    DOUBLE PRECISION :: DEL_X, DEL_Y    
    
    DEL_X = X2 - X1; DEL_Y = Y2 -Y1
    
    DO J = 0, M_NOW
        CALL AFFINE_MAPPING(GL_POINT_Y_T(J, PLY), &
                                        Y, Y1, DEL_Y)
                                        
        ! ON THE Y INTERFACE (BOTTOM)
        WRITE(5, FMT = '(F10.5, 2X, F10.5, 2X)', ADVANCE = 'NO') X1, Y 
        DO S = 1, NUM_OF_EQUATION
            WRITE(5, FMT = '(F10.5, 2X)', ADVANCE = 'NO') Q_Y1(J, S)
        ENDDO
        WRITE(5, FMT = '(3I6)') N_NOW, M_NOW, RANK                                
        
        ! INTERIOR NODES
        DO I = 0, N_NOW
            CALL AFFINE_MAPPING(GL_POINT_X_T(I, PLX), &
                                        X, X1, DEL_X)
                                        
            WRITE(5, FMT = '(F10.5, 2X, F10.5, 2X)', ADVANCE = 'NO') X, Y 
            DO S = 1, NUM_OF_EQUATION
                WRITE(5, FMT = '(F10.5, 2X)', ADVANCE = 'NO') Q(I, J, S)
            ENDDO
            WRITE(5, FMT = '(3I6)') N_NOW, M_NOW, RANK    
              
        ENDDO
        
        ! ON THE Y INTERFACE (TOP)
        WRITE(5, FMT = '(F10.5, 2X, F10.5, 2X)', ADVANCE = 'NO') X2, Y 
        DO S = 1, NUM_OF_EQUATION
            WRITE(5, FMT = '(F10.5, 2X)', ADVANCE = 'NO') Q_Y2(J, S)
        ENDDO
        WRITE(5, FMT = '(3I6)') N_NOW, M_NOW, RANK     
    
    ENDDO


END SUBROUTINE INTERIOR_MESH

END MODULE OUTPUT
