!-----------------------------------------------------------------------
!> @brief
!> Local data storage (inside each processor)
!-----------------------------------------------------------------------

MODULE LOCAL_STORAGE

USE MPI

IMPLICIT NONE 

INTEGER :: LOCAL_ELEM_NUM   !< LOCAL ELEMENT NUMBER

INTEGER :: ORIGINAL_ELEM_NUM    !< ELEMNT NUMBER BEFORE ADAPT

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: X_LOCAL   !< X COORDINATES
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: Y_LOCAL   !< Y COORDINATES



END MODULE LOCAL_STORAGE
