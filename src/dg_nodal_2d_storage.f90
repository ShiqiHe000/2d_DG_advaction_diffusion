!-----------------------------------------------------------------------
!> @brief
!> Storage for a nodal spectral method. Includes GL points, Gl weigths, 
!! first derivative matrices.
!-----------------------------------------------------------------------

MODULE NODAL_2D_STORAGE

USE MPI
USE PARAM, ONLY: N, M
USE BASIS

IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GL_POINT_X   !< GL POINTS IN X DIRECTION
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GL_POINT_Y   !< GL POINTS IN Y DIRECTION
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GL_W_X       !< GL WEIGHTS IN X DIRECTION
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GL_W_Y       !< GL WEIGHTS IN Y DIRECTION

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: FIRST_DER_X   !< FIRST DERIVATIVE IN X DIRECTION
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: FIRST_DER_Y   !< FIRST DERIVATIVE IN X DIRECTION

CONTAINS

SUBROUTINE GET_NODAL_2D_STORGAE

    IMPLICIT NONE
    
    ! ALLOCATE----------------------------------------------------------
    ALLOCATE(GL_POINT_X(0:N), GL_POINT_Y(0:M))
    ALLOCATE(GL_W_X(0:N), GL_W_Y(0:M))
    ALLOCATE(FIRST_DER_X(0:N, 0:N), FIRST_DER_Y(0:M, 0:M))
    !-------------------------------------------------------------------
    
    ! GET GL POINTS-----------------------------------------------------
    CALL GL(N, GL_POINT_X, GL_W_X)
    CALL GL(M, GL_POINT_Y, GL_W_Y)
    !-------------------------------------------------------------------

    ! GET FIRST DERIVATIVE MATRIX---------------------------------------
    CALL mth_Order_Polynomial_Derivative_Matrix(N, 1, GL_POINT_X, FIRST_DER_X)
    CALL mth_Order_Polynomial_Derivative_Matrix(M, 1, GL_POINT_Y, FIRST_DER_Y)
    !-------------------------------------------------------------------
    
!    PRINT *, GL_POINT_X
!    PRINT *, GL_W_Y


END SUBROUTINE GET_NODAL_2D_STORGAE


END MODULE NODAL_2D_STORAGE
