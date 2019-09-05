!-----------------------------------------------------------------------
!> @brief
!> Time derivative in 2d for DG approximation
!-----------------------------------------------------------------------

MODULE DG_TIME_DERIVATIVE

USE MPI
USE PARAM, ONLY: N, M, NUM_OF_EQUATION
USE RIEMANN_SOLVER
USE NODAL_2D_STORAGE
USE RIEMANN_SOLVER
USE FLUX_VECTORS

IMPLICIT NONE

CONTAINS

SUBROUTINE DG_TIME_DER 
!-----------------------------------------------------------------------
! Algorithm 93
!-----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: I, J, S
    
    DOUBLE PRECISION :: X, Y    ! AIMED POINT
    
    ALLOCATE(SOLUTION(0:N, 0:M, NUM_OF_EQUATION))
    
    ! X DIRECTION-------------------------------------------------------
    DO J=0, M
        Y = GL_POINT_Y(J)
        DO S=1, NUM_OF_EQUATION
            
        
        ENDDO
    
    ENDDO

END SUBROUTINE DG_TIME_DER


END MODULE DG_TIME_DERIVATIVE
