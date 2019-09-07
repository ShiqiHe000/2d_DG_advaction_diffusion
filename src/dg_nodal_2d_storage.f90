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
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: FIRST_DER_Y   !< FIRST DERIVATIVE IN Y DIRECTION


DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: M_FIRST_DER_X !< MODIFIED FIRST DERIVATIVE MATRIX
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: M_FIRST_DER_Y !< MODIFIED FIRST DERIVATIVE MATRIX

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: LAGRANGE_LEFT    !< LAGRANGE INTERPOLATING POLYNOMIAL ON THE LEFT BOUNDARY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: LAGRANGE_RIGHT   !< LAGRANGE INTERPOLATING POLYNOMIAL ON THE RIGHT BOUNDARY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: LAGRANGE_UP  !< LAGRANGE INTERPOLATING POLYNOMIAL ON THE UPPER BOUNDARY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: LAGRANGE_DOWN    !< LAGRANGE INTERPOLATING POLYNOMIAL ON THE LOWER BOUNDARY

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: NFLUX_X_L !< NUMERICAL FLUX AT LEFT BOUNDARY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: NFLUX_X_R !< NUMERICAL FLUX AT RIGHT BOUNDARY 
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: NFLUX_Y_D !< NUMERICAL FLUX AT BOTTOM  
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: NFLUX_Y_U !< NUMERICAL FLUX AT TOP 

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: FLUX_X !< FLUX 
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: FLUX_Y !< FLUX 

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: FLUX_DER_X    !< FLUX DERIVATIVE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: FLUX_DER_Y    !< FLUX DERIVATIVE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :, :) :: SOLUTION   !< SOLUTION(PRESSURE, VELOCITIES IN TWO DIRECTION)

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SOLUTION_INT_L   !< INTERIOR SOLUTION AT LEFT BOUNDARY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SOLUTION_INT_R   !< INTERIOR SOLUTION AT RIGHT BOUNDARY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SOLUTION_EXT_L   !< EXTERNAL SOLUTION AT LEFT BOUNDARY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SOLUTION_EXT_R   !< EXTERNAL SOLUTION AT RIGHT BOUNDARY

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :, :) :: SOLUTION_TIME_DER  !< SOLUTION TIME DERIVATIVE

CONTAINS

SUBROUTINE GET_NODAL_2D_STORGAE_BASIS
!-----------------------------------------------------------------------
! ALGORITHM 63
!> GL POINTS AND WEIGTHS, FIRST AND SECIND DERIVATIVE MATRICES
!-----------------------------------------------------------------------


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
!    print *, FIRST_DER_X


END SUBROUTINE GET_NODAL_2D_STORGAE_BASIS

SUBROUTINE GET_NODAL_2D_STORGAE_EXTENDS
!-----------------------------------------------------------------------
! ALGORITHM 89
!> LAGRANGE INTERPOLATES ON THE BOUNDARIES
!-----------------------------------------------------------------------

    IMPLICIT NONE
    
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BARY_X  ! BARYCENTRIC WEIGHTS
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BARY_Y  ! BARYCENTRIC WEIGHTS

    !-------------------------------------------------------------------
    ALLOCATE(LAGRANGE_LEFT(0:N), LAGRANGE_RIGHT(0:N))
    ALLOCATE(LAGRANGE_UP(0:M), LAGRANGE_DOWN(0:M))
    ALLOCATE(BARY_X(0:N), BARY_Y(0:M))
    !-------------------------------------------------------------------
    
    ! LEFT AND RIGHT BOUNDARY-------------------------------------------
    CALL BARW(N, GL_POINT_X, BARY_X)
    CALL LAGRANGE_INTERPOLATING_POLYNOMIAL(N, -1.0D0, GL_POINT_X, BARY_X, LAGRANGE_LEFT)
    CALL LAGRANGE_INTERPOLATING_POLYNOMIAL(N,  1.0D0, GL_POINT_X, BARY_X, LAGRANGE_RIGHT)
    !-------------------------------------------------------------------
    
    ! UPPER AND LOWER BOUNDARY------------------------------------------
    CALL BARW(M, GL_POINT_Y, BARY_Y)
    CALL LAGRANGE_INTERPOLATING_POLYNOMIAL(M, -1.0D0, GL_POINT_Y, BARY_Y, LAGRANGE_DOWN)
    CALL LAGRANGE_INTERPOLATING_POLYNOMIAL(M,  1.0D0, GL_POINT_Y, BARY_Y, LAGRANGE_UP)
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    DEALLOCATE(BARY_X, BARY_Y)
    !-------------------------------------------------------------------
    
!    PRINT *, LAGRANGE_LEFT
    
END SUBROUTINE GET_NODAL_2D_STORGAE_EXTENDS


END MODULE NODAL_2D_STORAGE
