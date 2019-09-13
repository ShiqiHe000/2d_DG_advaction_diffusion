!-----------------------------------------------------------------------
!> @brief
!> User defined initial conditions
!-----------------------------------------------------------------------

MODULE USER_DEFINES

USE MPI
USE PARAM, ONLY: K_X, K_Y, C

IMPLICIT NONE

DOUBLE PRECISION, PARAMETER :: D = 0.2D0 / (2.0D0 * DSQRT(DLOG(2.0D0)))   !< PARAMETER
DOUBLE PRECISION, PARAMETER :: X0 = -0.8D0  !< PARAMETER
DOUBLE PRECISION, PARAMETER :: Y0 = -0.8D0  !< PARAMETER

CONTAINS

SUBROUTINE INITIAL_CONDITION(N, M, N_EQUATIONS, Q, GL_X, GL_Y)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N, M !< POLY ORDER
    INTEGER, INTENT(IN) :: N_EQUATIONS  !< NUMBER OF EQUATION   
    
    INTEGER :: I, J
    
    DOUBLE PRECISION :: Q(0:N, 0:M, N_EQUATIONS)   !< INITIAL SOLUTION
    DOUBLE PRECISION :: GL_X(0:N), GL_Y(0:M)    !< GL_POINTS
    
    DOUBLE PRECISION :: INTER   !INTERMEDIATE
    DOUBLE PRECISION :: X, Y    ! COLLOCATION POINTS COORDINATES

    DO J=0, M
    
        Y = GL_Y(J)
    
        DO I=0, N
        
            X = GL_X(I)
            
            INTER = DEXP( - (K_X * (X - X0) + &
                            K_Y * (Y - Y0))**2 / D**2)
            
            Q(I, J, 1) = INTER 
            
            Q(I, J, 2) = K_X / C * INTER 
            
            Q(I, J, 3) = K_Y / C * INTER 
        
        ENDDO
    
    ENDDO

END SUBROUTINE INITIAL_CONDITION

SUBROUTINE EXACT_SOLUTION(N, M, N_EQUATIONS, GL_X, GL_Y, E, T)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N, M !< POLY ORDER
    INTEGER, INTENT(IN) :: N_EQUATIONS  !< NUMBER OF EQUATION 
    INTEGER :: I, J
    
    DOUBLE PRECISION :: E(0:N, 0:M, N_EQUATIONS)   !< EXACT SOLUTION AT TIME T
    DOUBLE PRECISION :: GL_X(0:N), GL_Y(0:M)    !< GL_POINTS
    DOUBLE PRECISION :: T   !< CURRENT TIME
    
    DOUBLE PRECISION :: INTER   !INTERMEDIATE
    DOUBLE PRECISION :: X, Y    ! COLLOCATION POINTS COORDINATES
    
    DO J=0, M
    
        Y = GL_Y(J)
    
        DO I=0, N
        
            X = GL_X(I)
            
            INTER = DEXP( - (K_X * (X - X0) + &
                            K_Y * (Y - Y0) - C * T)**2 / D**2)
            
            E(I, J, 1) = INTER 
            
            E(I, J, 2) = K_X / C * INTER 
            
            E(I, J, 3) = K_Y / C * INTER 
        
        ENDDO
    
    ENDDO


END SUBROUTINE EXACT_SOLUTION


SUBROUTINE INITIAL_SINUSOIDAL(N, M, N_EQUATIONS, Q, GL_X, GL_Y)

    USE PARAM, ONLY: PI

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N, M !< POLY ORDER
    INTEGER, INTENT(IN) :: N_EQUATIONS  !< NUMBER OF EQUATION   
    
    INTEGER :: I, J
    
    DOUBLE PRECISION :: Q(0:N, 0:M, N_EQUATIONS)   !< INITIAL SOLUTION
    DOUBLE PRECISION :: GL_X(0:N), GL_Y(0:M)    !< GL_POINTS
    
    DOUBLE PRECISION :: X, Y    ! COLLOCATION POINTS COORDINATES

    DO J=0, N
        X = GL_X(J)
        DO I=0, M
            Y = GL_Y(I)
            
            Q(I, J, 1) = -C * PI * DCOS(PI * X)
            Q(I, J, 2) = PI * DCOS(PI * X)
            Q(I, J, 3) = 0.0D0
        
        ENDDO
    
    ENDDO

END SUBROUTINE INITIAL_SINUSOIDAL

SUBROUTINE SIN_EXACT(N, M, N_EQUATIONS, GL_X, GL_Y, E, T)

    USE PARAM, ONLY: PI

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N, M !< POLY ORDER
    INTEGER, INTENT(IN) :: N_EQUATIONS  !< NUMBER OF EQUATION 
    INTEGER :: I, J
    
    DOUBLE PRECISION :: E(0:N, 0:M, N_EQUATIONS)   !< EXACT SOLUTION AT TIME T
    DOUBLE PRECISION :: GL_X(0:N), GL_Y(0:M)    !< GL_POINTS
    DOUBLE PRECISION :: T   !< CURRENT TIME
    
    DOUBLE PRECISION :: X, Y    ! COLLOCATION POINTS COORDINATES
    
    DO J=0, N
        X = GL_X(J)
        DO I=0, M
            Y = GL_Y(I)
            
            E(I, J, 1) = -C * PI * DCOS(PI * (X + C * T))
            E(I, J, 2) = PI * DCOS(PI * (X + C * T))
            E(I, J, 3) = 0.0D0
        
        ENDDO
    ENDDO


END SUBROUTINE SIN_EXACT

END MODULE USER_DEFINES
