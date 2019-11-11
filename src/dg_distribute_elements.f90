!-----------------------------------------------------------------------
!> @brief
!> After reading the mesh file, we spread the elements between processors
!! equally. We support element number < processors number, i.e., a 
!!processor could be assigned with zero elements.
!-----------------------------------------------------------------------

MODULE FIRST_DISTRIBUTE

USE MPI
USE NODAL_2D_STORAGE
USE PARAM
USE LOCAL_STORAGE

IMPLICIT NONE 

CONTAINS

!-----------------------------------------------------------------------
!> DISTRIBUTE ELEMENTS EQUALLY BETWEEN PROCESSORS.
!! SUPPORT PROCESSORS WITH ZERO ELEMENT.
!-----------------------------------------------------------------------
SUBROUTINE DISTRIBUTE_ELEM
    
    IMPLICIT NONE 
    
    INTEGER, DIMENSION(NUM_PROC) :: LOCAL_ELEM_NUMBER ! LOCAL
    INTEGER :: AVERAGE  ! AVERAGE LOAD
    INTEGER :: LAST     ! LAST PROCESSOR'S LOAD
    INTEGER :: I
    
    INTEGER, DIMENSION(NUM_PROC) :: SENDCOUNTS    !< Integer array (of length group size) specifying the number of elements to send to each processor.
    INTEGER, DIMENSION(NUM_PROC) :: DISPLS    !< Integer array (of length group size). Entry i specifies the displacement (relative to sendbuf) from which to take the outgoing data to process i. 
    
    LOCAL_ELEM_NUMBER = 0
    SENDCOUNTS = 0
    DISPLS = 0

    ALLOCATE(ELEM_RANGE(0:NUM_PROC))
    ELEM_RANGE = 0; ELEM_RANGE(0) = -1
        
    ! RANK0 COMPUTE THE AVERAGE LOAD------------------------------------
    IF (RANK == 0) THEN
        
        ORIGINAL_ELEM_NUM = NUM_OF_ELEMENT
        
        IF (NUM_OF_ELEMENT < NUM_PROC) THEN
        
            LOCAL_ELEM_NUMBER(1:NUM_OF_ELEMENT) = 1
            
            SENDCOUNTS(1:NUM_OF_ELEMENT) = 1
            
            ELEM_RANGE(1) = 0
            
        ELSE
            
            AVERAGE = NUM_OF_ELEMENT / NUM_PROC
            
            LAST = NUM_OF_ELEMENT - (NUM_PROC - 1) * AVERAGE
            
            ELEM_RANGE(1) = AVERAGE - 1 ! ELEMENT NUMBER START WITH 0 (HILBERT SCHEME) 
            
            LOCAL_ELEM_NUMBER(:) = AVERAGE
            
            SENDCOUNTS(:) = AVERAGE
            
            LOCAL_ELEM_NUMBER(NUM_PROC) = LAST
            
            SENDCOUNTS(NUM_PROC) = LAST
            
            
        ENDIF
        
        DO I = 2, NUM_PROC
            DISPLS(I) = LOCAL_ELEM_NUMBER(I-1) * 4     ! 4 COORDINATES
            ELEM_RANGE(I) = ELEM_RANGE(I-1) + LOCAL_ELEM_NUMBER(I)
        ENDDO
        
        
    ENDIF
    !-------------------------------------------------------------------
    
    ! SCATTER LOCAL ELEMENT NUMBER
    CALL MPI_SCATTER(LOCAL_ELEM_NUMBER, 1, MPI_INTEGER, LOCAL_ELEM_NUM, &
                    & 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERROR)
    
    ! BROADCAST TOTAL ELEMENT NUMBER
    CALL MPI_BCAST(ORIGINAL_ELEM_NUM, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERROR)
    
    ! ALLOCATE LOCAL STORAGE--------------------------------------------
    ALLOCATE(X_LOCAL(4, 0:LOCAL_ELEM_NUM-1))
    ALLOCATE(Y_LOCAL(4, 0:LOCAL_ELEM_NUM-1))
    
    IF(RANK /= 0) THEN
        ALLOCATE(X_HILBERT(4, ORIGINAL_ELEM_NUM))
        ALLOCATE(Y_HILBERT(4, ORIGINAL_ELEM_NUM))
    ENDIF
    !-------------------------------------------------------------------
    
    ! SCATTER DATA------------------------------------------------------
    CALL MPI_SCATTERV(X_HILBERT, SENDCOUNTS*4, DISPLS, MPI_DOUBLE_PRECISION, &
                        X_LOCAL, LOCAL_ELEM_NUM*4, MPI_DOUBLE_PRECISION, &
                        & 0, MPI_COMM_WORLD, IERROR)
                        
    CALL MPI_SCATTERV(Y_HILBERT, SENDCOUNTS*4, DISPLS, MPI_DOUBLE_PRECISION, &
                        Y_LOCAL, LOCAL_ELEM_NUM*4, MPI_DOUBLE_PRECISION, &
                        & 0, MPI_COMM_WORLD, IERROR)
                        
    CALL MPI_BCAST(ELEM_RANGE, NUM_PROC+1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERROR)
    ! ------------------------------------------------------------------
    
    DEALLOCATE(X_HILBERT, Y_HILBERT)
    IF (RANK == 0) THEN
        DEALLOCATE(DELTA_X, DELTA_Y)
    ENDIF
    
END SUBROUTINE DISTRIBUTE_ELEM


END MODULE FIRST_DISTRIBUTE
