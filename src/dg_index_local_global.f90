!-----------------------------------------------------------------------
!> @brief
!> Transform between local index and global index
!-----------------------------------------------------------------------

MODULE INDEX_LOCAL_GLOBAL

USE MPI
USE LOCAL_STORAGE

IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------
!> Input local rank and element global index, output local index.
!-----------------------------------------------------------------------
SUBROUTINE INDEX_GLOBAL_TO_LOCAL(RANK, G_INDEX, L_INDEX)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: RANK !< LOCAL RANK NUMBER
    INTEGER, INTENT(IN) :: G_INDEX  !< ELEMENT GLOBAL INDEX
    INTEGER, INTENT(OUT) :: L_INDEX  !< ELEMENT LOCAL INDEX
    
    L_INDEX = G_INDEX - (ELEM_RANGE(RANK) + 1)

END SUBROUTINE INDEX_GLOBAL_TO_LOCAL



END MODULE INDEX_LOCAL_GLOBAL



