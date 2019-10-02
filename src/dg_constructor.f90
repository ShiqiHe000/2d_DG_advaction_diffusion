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

IMPLICIT NONE

CONTAINS

SUBROUTINE CONSTRUCT_BASIS

    IMPLICIT NONE
   
   ! CONOSTRUCT DG BASIS-----------------------------------------------
    CALL CONSTRUCT_BASIS_STORAGE
   !-------------------------------------------------------------------
   
   ! INITIALIZE EACH ELEMENT WITH POLY LEVEL 1--------------------------
   ALLOCATE(PLEVEL_X(0:NUM_OF_ELEMENT-1))   ! element number start from 0 (use Hilbert scheme)
   ALLOCATE(PLEVEL_Y(0:NUM_OF_ELEMENT-1))
   
   PLEVEL_X = 1; PLEVEL_Y = 1
   !--------------------------------------------------------------------
    

END SUBROUTINE CONSTRUCT_BASIS


END MODULE DG_2D_CONSTRUCTOR
