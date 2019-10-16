!-----------------------------------------------------------------------
!> @brief
!> Map between an original set of point {x_j}_j ^N to 
!! a new set of points {X_j}_j ^M
!-----------------------------------------------------------------------

MODULE INTERPOLATE_TO_NEW_POINT

USE MPI
USE BASIS

IMPLICIT NONE 

CONTAINS

!-----------------------------------------------------------------------
! ALGORITHM 32
!> MATRIX FOR INTERPOLATION BETWEEN TWO SETS OF POINTS
!-----------------------------------------------------------------------
SUBROUTINE POLY_INTERPOLATION_MATRIX(N1, X, BW, N2, XI, T)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: N1   !< OLD POLYNOMIAL ORDER
    INTEGER, INTENT(IN) :: N2   !< NEW POLYNOMIAL ORDER
    
    INTEGER :: K, J
    
    DOUBLE PRECISION :: X(0:N1) !< ORIGINAL POINTS LOCATION
    DOUBLE PRECISION :: XI(0:N2)    !< NEW POINTS LOCATION
    DOUBLE PRECISION :: BW(0:N1)    !< BARYCENTRIC WEIGHTS
    
    DOUBLE PRECISION :: T(0:N2, 0:N1)   !< INTERPOLATE MATRIX (OUTPUT)
    
    DOUBLE PRECISION :: S, L   ! INTERMEDIATE VARIABLE
    
    LOGICAL :: ROW_MATCH
    LOGICAL :: FLAG
    
    DO K = 0, N2
        ROW_MATCH = .FALSE.
        
        DO J = 0, N1
        
            T(K, J) = 0.0D0
            
            CALL ALMOSTEQUAL(FLAG, XI(K), X(J))
            
            IF(FLAG) THEN
                ROW_MATCH = .TRUE.
                T(K, J) = 1.0D0
            ENDIF
        
        ENDDO
        
        IF(.NOT. FLAG) THEN
            S = 0.0D0
            
            DO J = 0, N1
                L = BW(J) / (XI(K) - X(J))
                T(K, J) = L
                S = S + L
            ENDDO
            
            DO J = 0, N1
                T(K, J) = T(K, J) / S
            ENDDO
            
        ENDIF
        
    
    ENDDO

END SUBROUTINE POLY_INTERPOLATION_MATRIX

!-----------------------------------------------------------------------
! Algorithm 33
!> Interpolation between two sets of points by matrix multipication.
!-----------------------------------------------------------------------
SUBROUTINE INTER_TO_NEW_POINT(N1, N2, T, F_OLD, F_NEW)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: N1   !< OLD POLY ORDER
    INTEGER, INTENT(IN) :: N2   !< NEW POLY ORDER
    
    INTEGER :: I, J
    
    DOUBLE PRECISION :: T(0:N2, 0:N1)   !< INTERPOLATE MATRIX
    DOUBLE PRECISION :: F_OLD(0:N1)     !< VALUE VECTOR OF OLD POINT SET
    DOUBLE PRECISION :: F_NEW(0:N2)     !< VALUE VECTOR OF NEW POINT SET(OUTPUT)
    
    DOUBLE PRECISION :: S
    
    DO I = 0, N2
        S = 0.0D0
        
        DO J = 0, N1
            S = S + T(I, J) * F_OLD(J)
        ENDDO
        
        F_NEW(I) = S
        
    ENDDO
    

END SUBROUTINE INTER_TO_NEW_POINT


SUBROUTINE INTERPOLATE_TO_NEW_SETS_POINTS(N1, N2, X, XI, M1, M2, Y, ETA, &
                                            F_OLD, F_NEW)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: N1   !< PLOY ORDER IN X DIRECTION (OLD)
    INTEGER, INTENT(IN) :: N2   !< PLOY ORDER IN X DIRECTION (NEW)
    INTEGER, INTENT(IN) :: M1   !< PLOY ORDER IN Y DIRECTION (OLD)
    INTEGER, INTENT(IN) :: M2   !< PLOY ORDER IN Y DIRECTION (NEW)
    
    INTEGER :: J, I
    
    DOUBLE PRECISION :: X(0:N1)     !< OLD NODES POSITION IN X DIRECTION
    DOUBLE PRECISION :: XI(0:N2)    !< NEW NODES POSITION IN X DIRECTION
    DOUBLE PRECISION :: Y(0:M1)     !< OLD NODES POSITION IN Y DIRECTION
    DOUBLE PRECISION :: ETA(0:M2)   !< NEW NODES POSITION IN Y DIRECTION
    
    DOUBLE PRECISION :: F_OLD(0:N1, 0:M1)   !< OLD VALUE ON THE OLD POINTS
    DOUBLE PRECISION :: F_NEW(0:N2, 0:M2)   !< NEW VALUE ON THE NEW POINTS
    DOUBLE PRECISION :: F_M(0:N2, 0:M1)     !< INTERMEDIATE MATRIX
    
    DOUBLE PRECISION :: BARYX(0:N1)  ! BARYCENTRIC WEIGHTS OF X
    DOUBLE PRECISION :: BARYY(0:M1)  ! BARYCENTRIC WEIGHTS OF Y
    
    DOUBLE PRECISION :: TX(0:N2, 0:N1)  ! INTERPOLATION MATRIX OF X
    DOUBLE PRECISION :: TY(0:M2, 0:M1)  ! INTERPOLATION MATRIX OF Y
    
    F_NEW = 0.0D0
    
    !-------------------------------------------------------------------
    CALL BARW(N1, X, BARYX)
    
    CALL POLY_INTERPOLATION_MATRIX(N1, X, BARYX, N2, XI, TX)
    
    DO J = 0, M1
        CALL INTER_TO_NEW_POINT(N1, N2, TX, F_OLD(:, J), F_M(:, J))
    ENDDO
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    CALL BARW(M1, Y, BARYY)
    
    CALL POLY_INTERPOLATION_MATRIX(M1, Y, BARYY, M2, ETA, TY)
    
    DO I = 0, N2
        CALL INTER_TO_NEW_POINT(M1, M2, TY, F_M(I, :), F_NEW(I, :))
    ENDDO
    !-------------------------------------------------------------------

END SUBROUTINE INTERPOLATE_TO_NEW_SETS_POINTS

END MODULE INTERPOLATE_TO_NEW_POINT
