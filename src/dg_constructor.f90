!-----------------------------------------------------------------------
!> @brief
!> A nodal DG 2D constructor.
!! This constructor computes the GL nodes and weights, Lagrange interpolating
!! polynomials at the boundaries, first order derivative matrices and 
!! modified them.
!-----------------------------------------------------------------------

MODULE DG_2D_CONSTRUCTOR

USE MPI
USE NODAL_2D_STORAGE
USE BASIS_STORAGE
USE LOCAL_STORAGE

IMPLICIT NONE

CONTAINS

SUBROUTINE CONSTRUCT_BASIS

    IMPLICIT NONE
   
   ! CONOSTRUCT DG BASIS-----------------------------------------------
    CALL CONSTRUCT_BASIS_STORAGE
   !-------------------------------------------------------------------
   
   
    ! INITIALIZE FUNDAMENTAL PARAMETERS---------------------------------
    ALLOCATE(A_LEVEL(0:LOCAL_ELEM_NUM-1))   ! element number start from 0 (use Hilbert scheme)
    A_LEVEL = 0

    ALLOCATE(PLEVEL_X(0:LOCAL_ELEM_NUM-1), PLEVEL_Y(0:LOCAL_ELEM_NUM-1))
    PLEVEL_X = 1; PLEVEL_Y = 1;
    !-------------------------------------------------------------------
    

END SUBROUTINE CONSTRUCT_BASIS


END MODULE DG_2D_CONSTRUCTOR
