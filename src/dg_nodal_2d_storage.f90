!-----------------------------------------------------------------------
!> @brief
!> Storage for a nodal spectral method. Includes GL points, Gl weigths, 
!! first derivative matrices.
!-----------------------------------------------------------------------

MODULE NODAL_2D_STORAGE

USE MPI
USE BASIS

IMPLICIT NONE

INTEGER :: NUM_OF_ELEMENT   !< TOTAL NUMBER OF ELEMENT

INTEGER :: NUM_OF_ELEMENT_X     !< NUMBER OF ELEMENT IN X DIRECTION
INTEGER :: NUM_OF_ELEMENT_Y     !< NUMBER OF ELEMENT IN Y DIRECTION

INTEGER, ALLOCATABLE, DIMENSION(:, :) :: DUAL_COORD   !< DUAL GRAPH COORDINATES

INTEGER, ALLOCATABLE, DIMENSION(:) :: PLEVEL_X     !< POLYNOMIAL LEVEL IN X(CAN BE TRANSFORM TO POLY ORDER)
INTEGER, ALLOCATABLE, DIMENSION(:) :: PLEVEL_Y     !< POLYNOMIAL LEVEL IN Y(CAN BE TRANSFORM TO POLY ORDER)

INTEGER :: LEVEL_MAX_X  !< MAXIMUM POLYNOMIAL LEVEL IN X
INTEGER :: LEVEL_MAX_Y  !< MAXIMUM POLYNOMIAL LEVEL IN Y

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: ELEM_X_POSITION    !< ELEMENT K'S X COORDINATES IN DESIGNED SEQUENCE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: ELEM_Y_POSITION    !< ELEMENT K'S Y COORDINATES IN DESIGNED SEQUENCE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :)  :: X_HILBERT    !< X COORDINATE OF ELEMENT K ON HILBER CURVE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :)  :: Y_HILBERT    !< Y COORDINATE OF ELEMENT K ON HILBER CURVE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DELTA_X !< ELEMENT SIZE IN X DIRECTION
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DELTA_Y !< ELEMENT SIZE IN Y DIRECTION

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: GL_POINT_X_T   !< GL POINTS TABLE IN X DIRECTION
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: GL_POINT_Y_T   !< GL POINTS TABLE IN Y DIRECTION
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: GL_W_X_T       !< GL WEIGHTS TABLE IN X DIRECTION
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: GL_W_Y_T       !< GL WEIGHTS TABLE IN Y DIRECTION

!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GL_POINT_X   !< GL POINTS IN X DIRECTION
!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GL_POINT_Y   !< GL POINTS IN Y DIRECTION
!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GL_W_X       !< GL WEIGHTS IN X DIRECTION
!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GL_W_Y       !< GL WEIGHTS IN Y DIRECTION

!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: FIRST_DER_X   !< FIRST DERIVATIVE IN X DIRECTION
!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: FIRST_DER_Y   !< FIRST DERIVATIVE IN Y DIRECTION

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :, :) :: FIRST_DER_X_T   !< FIRST DERIVATIVE TABLE IN X DIRECTION
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :, :) :: FIRST_DER_Y_T   !< FIRST DERIVATIVE TABLE IN Y DIRECTION


!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: M_FIRST_DER_X !< MODIFIED FIRST DERIVATIVE MATRIX
!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: M_FIRST_DER_Y !< MODIFIED FIRST DERIVATIVE MATRIX

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :, :) :: M_FIRST_DER_X_T !< MODIFIED FIRST DERIVATIVE MATRIX TABLE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :, :) :: M_FIRST_DER_Y_T !< MODIFIED FIRST DERIVATIVE MATRIX TABLE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: LAGRANGE_LEFT_T   !< LAGRANGE INTERPOLATING POLYNOMIAL TABLE ON THE LEFT BOUNDARY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: LAGRANGE_RIGHT_T  !< LAGRANGE INTERPOLATING POLYNOMIAL TABLE ON THE RIGHT BOUNDARY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: LAGRANGE_UP_T  !< LAGRANGE INTERPOLATING POLYNOMIAL TABLE ON THE UPPER BOUNDARY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: LAGRANGE_DOWN_T    !< LAGRANGE INTERPOLATING POLYNOMIAL TABLE ON THE LOWER BOUNDARY

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: NFLUX_X_L !< NUMERICAL FLUX AT LEFT BOUNDARY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: NFLUX_X_R !< NUMERICAL FLUX AT RIGHT BOUNDARY 
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: NFLUX_Y_D !< NUMERICAL FLUX AT BOTTOM  
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: NFLUX_Y_U !< NUMERICAL FLUX AT TOP 

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: FLUX_X !< FLUX 
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: FLUX_Y !< FLUX 

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: FLUX_DER_X    !< FLUX DERIVATIVE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: FLUX_DER_Y    !< FLUX DERIVATIVE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :, :, :) :: SOLUTION   !< SOLUTION(PRESSURE, VELOCITIES IN TWO DIRECTION)

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SOLUTION_INT_L   !< INTERIOR SOLUTION AT LEFT BOUNDARY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SOLUTION_INT_R   !< INTERIOR SOLUTION AT RIGHT BOUNDARY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SOLUTION_EXT_L   !< EXTERNAL SOLUTION AT LEFT BOUNDARY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SOLUTION_EXT_R   !< EXTERNAL SOLUTION AT RIGHT BOUNDARY

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :, :) :: SOLUTION_TIME_DER  !< SOLUTION TIME DERIVATIVE

CONTAINS

!-----------------------------------------------------------------------
! ALGORITHM 63
!> GL POINTS AND WEIGTHS, FIRST AND SECIND DERIVATIVE MATRICES
!-----------------------------------------------------------------------
SUBROUTINE GET_NODAL_2D_STORAGE_BASIS(N1, GL_P, GL_WEIGHT, &
                                      FIRST_DER)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N1   !< POLY ORDER 
    
    DOUBLE PRECISION :: GL_P(0:N1)        !< GL POINTS 
    DOUBLE PRECISION :: GL_WEIGHT(0:N1)   !< GL WEIGHTS
    DOUBLE PRECISION :: FIRST_DER(0:N1, 0:N1)   ! FIRST DIREVATIVE MATRIX
    
    ! GET GL POINTS-----------------------------------------------------
    CALL GL(N1, GL_P, GL_WEIGHT)
    !-------------------------------------------------------------------

    ! GET FIRST DERIVATIVE MATRIX--------------------------------------- 
    CALL mth_Order_Polynomial_Derivative_Matrix(N1, 1, GL_P, FIRST_DER)
    !-------------------------------------------------------------------


END SUBROUTINE GET_NODAL_2D_STORAGE_BASIS

!-----------------------------------------------------------------------
! ALGORITHM 89
!> LAGRANGE INTERPOLATES ON THE BOUNDARIES
!-----------------------------------------------------------------------
SUBROUTINE GET_NODAL_2D_STORAGE_EXTENDS(N1, LAG_L, LAG_R, GL_P)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N1   !< POLY ORDER
    
    DOUBLE PRECISION :: LAG_L(0:N1) !< LAGRANGE INTERPOLATE ON THE LEFT BOUNDARY
    DOUBLE PRECISION :: LAG_R(0:N1) !< LAGRANGE INTERPOLATE ON THE RIGHT BOUNDARY
    
    DOUBLE PRECISION :: GL_P(0:N1)        !< GL POINTS 
    
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BARY_X  ! BARYCENTRIC WEIGHTS

    !-------------------------------------------------------------------
    ALLOCATE(BARY_X(0:N1))
    !-------------------------------------------------------------------
    
    ! LEFT AND RIGHT BOUNDARY-------------------------------------------
    CALL BARW(N1, GL_P, BARY_X)
    CALL LAGRANGE_INTERPOLATING_POLYNOMIAL(N1, -1.0D0, GL_P, BARY_X, LAG_L)
    CALL LAGRANGE_INTERPOLATING_POLYNOMIAL(N1,  1.0D0, GL_P, BARY_X, LAG_R)
    !-------------------------------------------------------------------
    
    ! UPPER AND LOWER BOUNDARY------------------------------------------
!    CALL BARW(M, GL_POINT_Y, BARY_Y)
!    CALL LAGRANGE_INTERPOLATING_POLYNOMIAL(M, -1.0D0, GL_POINT_Y, BARY_Y, LAGRANGE_DOWN)
!    CALL LAGRANGE_INTERPOLATING_POLYNOMIAL(M,  1.0D0, GL_POINT_Y, BARY_Y, LAGRANGE_UP)
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    DEALLOCATE(BARY_X)
    !-------------------------------------------------------------------

    
END SUBROUTINE GET_NODAL_2D_STORAGE_EXTENDS


END MODULE NODAL_2D_STORAGE
