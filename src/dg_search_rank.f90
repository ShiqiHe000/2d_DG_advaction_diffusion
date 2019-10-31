!-----------------------------------------------------------------------
!> @brief 
!> This module search the corresponding rank of the input element number
!-----------------------------------------------------------------------

MODULE SEARCH_RANK

USE MPI
USE LOCAL_STORAGE
USE PARAM

IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------
!> Input the target number (root), returns the rank of processor it belongs
!-----------------------------------------------------------------------
SUBROUTINE FIND_RANK(TARGET_NUM, TARGET_RANK)

    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: TARGET_NUM   !< TARGET ELEMENT NUMBER
    INTEGER :: TARGET_RANK  !< RANK OF TARGET
    INTEGER :: S
    
    TARGET_RANK = RANK

    IF(TARGET_NUM > ELEM_RANGE(RANK + 1)) THEN
    
        TARGET_RANK = TARGET_RANK + 1
        
        DO S = RANK+2, NUM_PROC
            IF (TARGET_NUM > ELEM_RANGE(S)) THEN
                TARGET_RANK = TARGET_RANK + 1
            ELSE
                RETURN
            ENDIF
        ENDDO
        
    ELSE
        TARGET_RANK = TARGET_RANK - 1
    
        DO S = RANK-1, 1, -1
            IF (TARGET_NUM <= ELEM_RANGE(S)) THEN
                TARGET_RANK = TARGET_RANK -1
            ELSE
                RETURN
            ENDIF
        ENDDO
    ENDIF
    

END SUBROUTINE FIND_RANK

END MODULE SEARCH_RANK
