!-----------------------------------------------------------------------
!> @brief
!> A nodal DG 2D constructor.
!! This constructor computes the GL nodes and weights, Lagrange interpolating
!! polynomials at the boundaries, first order derivative matrices and 
!! modified them.
!-----------------------------------------------------------------------

MODULE DG_2D_CONSTRUCTOR

USE MPI
USE PARAM, ONLY: N, M, X_L, X_R, Y_L, Y_R
USE NODAL_2D_STORAGE

IMPLICIT NONE

CONTAINS

SUBROUTINE CONSTRUCT_BASIS
!-----------------------------------------------------------------------
! Algorithm 91
!-----------------------------------------------------------------------

    IMPLICIT NONE
    
    INTEGER :: I, J
    
    !-------------------------------------------------------------------
    ALLOCATE(M_FIRST_DER_X(0:N, 0:N), M_FIRST_DER_Y(0:M, 0:M))
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    CALL GET_NODAL_2D_STORGAE_BASIS
    
    CALL GET_NODAL_2D_STORGAE_EXTENDS
    !-------------------------------------------------------------------
    
    ! MODIFY THE FIRST ORDER DERIVATIVE ATRICES-------------------------
    DO J=0, N
        DO I=0, N
            ! X DIRECTION
            M_FIRST_DER_X(I, J) = - FIRST_DER_X(J, I) * GL_W_X(J) / GL_W_X(I)
        ENDDO
    ENDDO
    
    DO J=0, M
        DO I=0, M
            ! Y DIRECTION
            M_FIRST_DER_Y(I, J) = - FIRST_DER_Y(J, I) * GL_W_Y(J) / GL_W_Y(I)
        ENDDO
    ENDDO
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    DEALLOCATE(FIRST_DER_X, FIRST_DER_Y)
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    DELTA_X = X_R - X_L
    DELTA_Y = Y_R - Y_L
    !-------------------------------------------------------------------

END SUBROUTINE CONSTRUCT_BASIS


END MODULE DG_2D_CONSTRUCTOR
