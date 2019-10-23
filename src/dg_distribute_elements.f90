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
    
    INTEGER, ALLOCATABLE, DIMENSION(:) :: LOCAL_ELEM_NUMBER ! LOCAL
    INTEGER :: AVERAGE  ! AVERAGE LOAD
    INTEGER :: LAST     ! LAST PROCESSOR'S LOAD
    INTEGER :: I
    
    INTEGER, ALLOCATABLE, DIMENSION(:) :: SENDCOUNTS    !< Integer array (of length group size) specifying the number of elements to send to each processor.
    INTEGER, ALLOCATABLE, DIMENSION(:) :: DISPLS    !< Integer array (of length group size). Entry i specifies the displacement (relative to sendbuf) from which to take the outgoing data to process i. 

    ! RANK0 COMPUTE THE AVERAGE LOAD------------------------------------
    IF (RANK == 0) THEN
    
        ALLOCATE(SENDCOUNTS(NUM_PROC))
        SENDCOUNTS = 0
        
        ALLOCATE(DISPLS(NUM_PROC))
        DISPLS = 0
    
        ORIGINAL_ELEM_NUM = NUM_OF_ELEMENT
        
        IF (NUM_OF_ELEMENT < NUM_PROC) THEN
        
            ALLOCATE(LOCAL_ELEM_NUMBER(NUM_PROC))
            LOCAL_ELEM_NUMBER = 0
            
            LOCAL_ELEM_NUMBER(1:NUM_OF_ELEMENT) = 1
            
            SENDCOUNTS(1:NUM_OF_ELEMENT) = 1
            
        ELSE
            
            ALLOCATE(LOCAL_ELEM_NUMBER(NUM_PROC))
            LOCAL_ELEM_NUMBER = 0
            
            AVERAGE = NUM_OF_ELEMENT / NUM_PROC
            
            LAST = NUM_OF_ELEMENT - (NUM_PROC - 1) * AVERAGE
            
            LOCAL_ELEM_NUMBER(:) = AVERAGE
            
            SENDCOUNTS(:) = AVERAGE
            
            LOCAL_ELEM_NUMBER(NUM_PROC) = LAST
            
            SENDCOUNTS(NUM_PROC) = LAST
            
            
        ENDIF
        
        DO I = 2, NUM_PROC
            DISPLS(I) = LOCAL_ELEM_NUMBER(I-1) * 4     ! 4 COORDINATES
        ENDDO
        
    ENDIF
    !-------------------------------------------------------------------
    
    ! SCATTER LOCAL ELEMENT NUMBER
    CALL MPI_SCATTER(LOCAL_ELEM_NUMBER, 1, MPI_INT, LOCAL_ELEM_NUM, &
                    & 1, MPI_INT, 0, MPI_COMM_WORLD, IERROR)
    
    ! BROADCAST TOTAL ELEMENT NUMBER
    CALL MPI_BCAST(ORIGINAL_ELEM_NUM, 1, MPI_INT, 0, MPI_COMM_WORLD, IERROR)
    
    ! ALLOCATE LOCAL STORAGE--------------------------------------------
    ALLOCATE(X_LOCAL(4, LOCAL_ELEM_NUM))
    ALLOCATE(Y_LOCAL(4, LOCAL_ELEM_NUM))
    !-------------------------------------------------------------------
    
    ! SCATTER DATA------------------------------------------------------
    CALL MPI_SCATTERV(X_HILBERT, SENDCOUNTS*4, DISPLS, MPI_DOUBLE_PRECISION, &
                        X_LOCAL, LOCAL_ELEM_NUM*4, MPI_DOUBLE_PRECISION, &
                        & 0, MPI_COMM_WORLD, IERROR)
                        
    CALL MPI_SCATTERV(Y_HILBERT, SENDCOUNTS*4, DISPLS, MPI_DOUBLE_PRECISION, &
                        Y_LOCAL, LOCAL_ELEM_NUM*4, MPI_DOUBLE_PRECISION, &
                        & 0, MPI_COMM_WORLD, IERROR)
    ! ------------------------------------------------------------------

END SUBROUTINE DISTRIBUTE_ELEM


END MODULE FIRST_DISTRIBUTE
