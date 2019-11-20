!-----------------------------------------------------------------------
!> @brief
!> Exchage ghost layger information after we construct the element 
!! interfaces.
!-----------------------------------------------------------------------

MODULE MESSAGE_EXCHAGE

USE MPI
USE LOCAL_STORAGE
USE INDEX_LOCAL_GLOBAL
USE PARAM
USE hilbert
USE SEARCH_RANK

IMPLICIT NONE 

CONTAINS

!-----------------------------------------------------------------------
!> EXCHANGE SOLUTIONS ON MPI INTERFACES (ALL ELEMENTS).
!! SYNCHRONIZE IN THE END.
!-----------------------------------------------------------------------
SUBROUTINE EXCHANGE_SOLUTION_X
    
    IMPLICIT NONE 
    
    INTEGER, DIMENSION(LOCAL_ELEM_NUM) :: ARRAY_OF_REQUESTS_SEND ! Array of requests (array of handles). (INPUT)
    INTEGER, DIMENSION(MPI_STATUS_SIZE, LOCAL_ELEM_NUM) :: ARRAY_OF_STATUSES_S ! Array of status objects (array of status). (OUTPUT)
    
    INTEGER :: NUM_ISEND
    INTEGER :: IERROR
    INTEGER :: K
    
    ARRAY_OF_REQUESTS_SEND = 0
    NUM_ISEND = 0
    
    DO K = 0, LOCAL_ELEM_NUM-1
        CALL NON_BLOCKING_EXCAHNGE_INTERFACE_X(MPI_B_FLAG(1, K), &
                                               MPI_B_FLAG(2, K), K, &
                                               NUM_ISEND, &
                                               ARRAY_OF_REQUESTS_SEND, &
                                               SOLUTION_INT_R(:, :, K), &
                                               GHOST(:, :, K))
        
    ENDDO
    
    CALL MPI_WAITALL(NUM_ISEND, ARRAY_OF_REQUESTS_SEND(1:NUM_ISEND), &
                    ARRAY_OF_STATUSES_S(:, 1:NUM_ISEND), IERROR)

END SUBROUTINE EXCHANGE_SOLUTION_X

!-----------------------------------------------------------------------
!> EXCHANGE NUMERICAL FLUXES ON MPI INTERFACES (ALL ELEMENTS).
!! SYNCHRONIZE IN THE END. 
!-----------------------------------------------------------------------
SUBROUTINE EXCHANGE_NFLUX_X
    
    IMPLICIT NONE 
    
    INTEGER, DIMENSION(LOCAL_ELEM_NUM) :: ARRAY_OF_REQUESTS_SEND ! Array of requests (array of handles). (INPUT)
    INTEGER, DIMENSION(MPI_STATUS_SIZE, LOCAL_ELEM_NUM) :: ARRAY_OF_STATUSES_S ! Array of status objects (array of status). (OUTPUT)
    
    INTEGER :: NUM_ISEND
    INTEGER :: IERROR
    INTEGER :: K
    
    ARRAY_OF_REQUESTS_SEND = 0
    NUM_ISEND = 0
    
    DO K = 0, LOCAL_ELEM_NUM-1
        CALL NON_BLOCKING_EXCAHNGE_NFLUX_X(MPI_B_FLAG(1, K), &
                                               MPI_B_FLAG(2, K), K, &
                                               NUM_ISEND, &
                                               ARRAY_OF_REQUESTS_SEND, &
                                               -NFLUX_X_L(:, :, K), &
                                               NFLUX_X_R(:, :, K))
        
    ENDDO
    
    CALL MPI_WAITALL(NUM_ISEND, ARRAY_OF_REQUESTS_SEND(1:NUM_ISEND), &
                    ARRAY_OF_STATUSES_S(:, 1:NUM_ISEND), IERROR)

END SUBROUTINE EXCHANGE_NFLUX_X

!-----------------------------------------------------------------------
!> INSIDE ONE ELEMENT USE NON-BLOCKING COMMUNICATION TO UPDATE MPI INTERFACES
!! X DIRECTION
!-----------------------------------------------------------------------
SUBROUTINE NON_BLOCKING_EXCAHNGE_INTERFACE_X(FLAG1, FLAG2, LELEM_K, &
                                            NUM_SEND, &
                                            ARRAY_OF_REQUESTS_SEND, &
                                            SEND_BUFF, RECV_BUFF)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: LELEM_K  !< LOCAL ELEMENT NUMBER (START WITH 0)
    INTEGER :: RELEM_K  ! ROOT ELEMENT NUMBER (START WITH 0)
    INTEGER :: I, J     ! ELEMENT COORDINATES
    INTEGER :: NEIGHBOUR_R  ! NEIGHBOUR ELEMENT NUMBER (ROOT)
    INTEGER :: TARGET_RANK  ! NEIGHBOUR'S RANK
    INTEGER :: SR_COUNT   ! Number of elements in send/recv buffer 
    INTEGER :: REQUEST  ! Communication request (handle).  
    INTEGER :: IERROR
    INTEGER :: NUM_SEND ! COUND ISEND TIMES
    INTEGER :: ARRAY_OF_REQUESTS_SEND(LOCAL_ELEM_NUM) ! RECORD COMMINICATION REQUEST HANDLE
    INTEGER :: STATUS(MPI_STATUS_SIZE)

    LOGICAL, INTENT(IN) :: FLAG1 !< MPI BOUNDARY FLAG (FACE1/3)
    LOGICAL, INTENT(IN) :: FLAG2 !< MPI BOUNDARY FLAG (FACE2/4)
    
    DOUBLE PRECISION :: SEND_BUFF(0:MMAX, NUM_OF_EQUATION)  !< SEND BUFFER
    DOUBLE PRECISION :: RECV_BUFF(0:MMAX, NUM_OF_EQUATION)  !< RECV BUFFER
    
    SR_COUNT = (1 + MMAX) * NUM_OF_EQUATION

    CALL INDEX_LOCAL_TO_GLOBAL(RANK, LELEM_K, RELEM_K)  ! CONVERT ELEMENT LOCAL INDEX TO GLOBAL INDEX

    IF(FLAG1) THEN  ! RECV INFORMATION REMOTELY
        
        ! USE MPI WILDCARDS: MPI_ANY_TAG, MPI_ANY_SOURCE
        CALL MPI_RECV(RECV_BUFF, SR_COUNT, MPI_DOUBLE_PRECISION, &
                        MPI_ANY_SOURCE, RELEM_K, &
                        MPI_COMM_WORLD, STATUS, IERROR)
        
    ENDIF
    
    IF(FLAG2) THEN
    
        CALL d2xy ( EXP_X, RELEM_K, J, I )  ! ELEMENT COORDINATE
        
        CALL xy2d ( EXP_X, J, I+1, NEIGHBOUR_R )    ! NEIGHBOUR ROOT ELEMENT NUMBER
        
        CALL FIND_RANK(NEIGHBOUR_R, TARGET_RANK)     ! ENIGHBOUR'S RANK
    
        CALL MPI_ISEND(SEND_BUFF, SR_COUNT, MPI_DOUBLE_PRECISION, &
                        TARGET_RANK, NEIGHBOUR_R, MPI_COMM_WORLD, &
                        REQUEST, IERROR)
        
        NUM_SEND = NUM_SEND + 1
        
        ARRAY_OF_REQUESTS_SEND(NUM_SEND) = REQUEST
        
    ENDIF
    
END SUBROUTINE NON_BLOCKING_EXCAHNGE_INTERFACE_X


!-----------------------------------------------------------------------
!> INSIDE ONE ELEMENT USE NON-BLOCKING COMMUNICATION TO UPDATE MPI INTERFACES
!! X DIRECTION
!-----------------------------------------------------------------------
SUBROUTINE NON_BLOCKING_EXCAHNGE_NFLUX_X(FLAG1, FLAG2, LELEM_K, &
                                            NUM_SEND, &
                                            ARRAY_OF_REQUESTS_SEND, &
                                            SEND_BUFF, RECV_BUFF)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: LELEM_K  !< LOCAL ELEMENT NUMBER (START WITH 0)
    INTEGER :: RELEM_K  ! ROOT ELEMENT NUMBER (START WITH 0)
    INTEGER :: I, J     ! ELEMENT COORDINATES
    INTEGER :: NEIGHBOUR_R  ! NEIGHBOUR ELEMENT NUMBER (ROOT)
    INTEGER :: TARGET_RANK  ! NEIGHBOUR'S RANK
    INTEGER :: SR_COUNT   ! Number of elements in send/recv buffer 
    INTEGER :: REQUEST  ! Communication request (handle).  
    INTEGER :: IERROR
    INTEGER :: NUM_SEND ! COUND ISEND TIMES
    INTEGER :: ARRAY_OF_REQUESTS_SEND(LOCAL_ELEM_NUM) ! RECORD COMMINICATION REQUEST HANDLE
    INTEGER :: STATUS(MPI_STATUS_SIZE)

    LOGICAL, INTENT(IN) :: FLAG1 !< MPI BOUNDARY FLAG (FACE1/3)
    LOGICAL, INTENT(IN) :: FLAG2 !< MPI BOUNDARY FLAG (FACE2/4)
    
    DOUBLE PRECISION :: SEND_BUFF(0:NMAX, NUM_OF_EQUATION)  !< SEND BUFFER
    DOUBLE PRECISION :: RECV_BUFF(0:NMAX, NUM_OF_EQUATION)  !< RECV BUFFER
    
    SR_COUNT = (1 + NMAX) * NUM_OF_EQUATION
    
    CALL INDEX_LOCAL_TO_GLOBAL(RANK, LELEM_K, RELEM_K)  ! CONVERT ELEMENT LOCAL INDEX TO GLOBAL INDEX

    IF(FLAG2) THEN  ! RECV INFORMATION REMOTELY
    
        ! USE MPI WILDCARDS: MPI_ANY_TAG, MPI_ANY_SOURCE
        CALL MPI_RECV(RECV_BUFF, SR_COUNT, MPI_DOUBLE_PRECISION, &
                        MPI_ANY_SOURCE, RELEM_K, &
                        MPI_COMM_WORLD, STATUS, IERROR)
        
    ENDIF
    
    IF(FLAG1) THEN
    
        CALL d2xy ( EXP_X, RELEM_K, J, I )  ! ELEMENT COORDINATE
        
        CALL xy2d ( EXP_X, J, I-1, NEIGHBOUR_R )    ! NEIGHBOUR ROOT ELEMENT NUMBER
        
        CALL FIND_RANK(NEIGHBOUR_R, TARGET_RANK)     ! ENIGHBOUR'S RANK
    
        CALL MPI_ISEND(SEND_BUFF, SR_COUNT, MPI_DOUBLE_PRECISION, &
                        TARGET_RANK, NEIGHBOUR_R, MPI_COMM_WORLD, &
                        REQUEST, IERROR)
        
        NUM_SEND = NUM_SEND + 1
        
        ARRAY_OF_REQUESTS_SEND(NUM_SEND) = REQUEST
        
    ENDIF


END SUBROUTINE NON_BLOCKING_EXCAHNGE_NFLUX_X


!-----------------------------------------------------------------------
!> EXCHANGE SOLUTIONS ON MPI INTERFACES (ALL ELEMENTS).
!! SYNCHRONIZE IN THE END.
!-----------------------------------------------------------------------
SUBROUTINE EXCHANGE_SOLUTION_Y
    
    IMPLICIT NONE 
    
    INTEGER, DIMENSION(LOCAL_ELEM_NUM) :: ARRAY_OF_REQUESTS_SEND ! Array of requests (array of handles). (INPUT)
    INTEGER, DIMENSION(MPI_STATUS_SIZE, LOCAL_ELEM_NUM) :: ARRAY_OF_STATUSES_S ! Array of status objects (array of status). (OUTPUT)
    
    INTEGER :: NUM_ISEND
    INTEGER :: IERROR
    INTEGER :: K
    
    ARRAY_OF_REQUESTS_SEND = 0
    NUM_ISEND = 0
    
    DO K = 0, LOCAL_ELEM_NUM-1
        CALL NON_BLOCKING_EXCAHNGE_INTERFACE_Y(MPI_B_FLAG(3, K), &
                                               MPI_B_FLAG(4, K), K, &
                                               NUM_ISEND, &
                                               ARRAY_OF_REQUESTS_SEND, &
                                               SOLUTION_INT_R(:, :, K), &
                                               GHOST(:, :, K))
        
    ENDDO
    
    CALL MPI_WAITALL(NUM_ISEND, ARRAY_OF_REQUESTS_SEND(1:NUM_ISEND), &
                    ARRAY_OF_STATUSES_S(:, 1:NUM_ISEND), IERROR)

END SUBROUTINE EXCHANGE_SOLUTION_Y

!-----------------------------------------------------------------------
!> EXCHANGE NUMERICAL FLUXES ON MPI INTERFACES (ALL ELEMENTS).
!! SYNCHRONIZE IN THE END. 
!-----------------------------------------------------------------------
SUBROUTINE EXCHANGE_NFLUX_Y
    
    IMPLICIT NONE 
    
    INTEGER, DIMENSION(LOCAL_ELEM_NUM) :: ARRAY_OF_REQUESTS_SEND ! Array of requests (array of handles). (INPUT)
    INTEGER, DIMENSION(MPI_STATUS_SIZE, LOCAL_ELEM_NUM) :: ARRAY_OF_STATUSES_S ! Array of status objects (array of status). (OUTPUT)
    
    INTEGER :: NUM_ISEND
    INTEGER :: IERROR
    INTEGER :: K
    
    ARRAY_OF_REQUESTS_SEND = 0
    NUM_ISEND = 0
    
    DO K = 0, LOCAL_ELEM_NUM-1
        CALL NON_BLOCKING_EXCAHNGE_NFLUX_Y(MPI_B_FLAG(3, K), &
                                               MPI_B_FLAG(4, K), K, &
                                               NUM_ISEND, &
                                               ARRAY_OF_REQUESTS_SEND, &
                                               -NFLUX_Y_D(:, :, K), &
                                               NFLUX_Y_U(:, :, K))
        
    ENDDO
    
    CALL MPI_WAITALL(NUM_ISEND, ARRAY_OF_REQUESTS_SEND(1:NUM_ISEND), &
                    ARRAY_OF_STATUSES_S(:, 1:NUM_ISEND), IERROR)

END SUBROUTINE EXCHANGE_NFLUX_Y


!-----------------------------------------------------------------------
!> INSIDE ONE ELEMENT USE NON-BLOCKING COMMUNICATION TO UPDATE MPI INTERFACES
!! Y DIRECTION
!-----------------------------------------------------------------------
SUBROUTINE NON_BLOCKING_EXCAHNGE_INTERFACE_Y(FLAG1, FLAG2, LELEM_K, &
                                            NUM_SEND, &
                                            ARRAY_OF_REQUESTS_SEND, &
                                            SEND_BUFF, RECV_BUFF)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: LELEM_K  !< LOCAL ELEMENT NUMBER (START WITH 0)
    INTEGER :: RELEM_K  ! ROOT ELEMENT NUMBER (START WITH 0)
    INTEGER :: I, J     ! ELEMENT COORDINATES
    INTEGER :: NEIGHBOUR_R  ! NEIGHBOUR ELEMENT NUMBER (ROOT)
    INTEGER :: TARGET_RANK  ! NEIGHBOUR'S RANK
    INTEGER :: SR_COUNT   ! Number of elements in send/recv buffer 
    INTEGER :: REQUEST  ! Communication request (handle).  
    INTEGER :: IERROR
    INTEGER :: NUM_SEND ! COUND ISEND TIMES
    INTEGER :: ARRAY_OF_REQUESTS_SEND(LOCAL_ELEM_NUM) ! RECORD COMMINICATION REQUEST HANDLE
    INTEGER :: STATUS(MPI_STATUS_SIZE)

    LOGICAL, INTENT(IN) :: FLAG1 !< MPI BOUNDARY FLAG (FACE1/3)
    LOGICAL, INTENT(IN) :: FLAG2 !< MPI BOUNDARY FLAG (FACE2/4)
    
    DOUBLE PRECISION :: SEND_BUFF(0:NMAX, NUM_OF_EQUATION)  !< SEND BUFFER
    DOUBLE PRECISION :: RECV_BUFF(0:NMAX, NUM_OF_EQUATION)  !< SEND BUFFER
    
    SR_COUNT = (1 + NMAX) * NUM_OF_EQUATION
    
    CALL INDEX_LOCAL_TO_GLOBAL(RANK, LELEM_K, RELEM_K)  ! CONVERT ELEMENT LOCAL INDEX TO GLOBAL INDEX

    IF(FLAG1) THEN  ! RECV INFORMATION REMOTELY
    
        ! USE MPI WILDCARDS: MPI_ANY_TAG, MPI_ANY_SOURCE
        CALL MPI_RECV(RECV_BUFF, SR_COUNT, MPI_DOUBLE_PRECISION, &
                        MPI_ANY_SOURCE, RELEM_K, &
                        MPI_COMM_WORLD, STATUS, IERROR)
        
    ENDIF
    
    IF(FLAG2) THEN
    
        CALL d2xy ( EXP_X, RELEM_K, J, I )  ! ELEMENT COORDINATE
        
        CALL xy2d ( EXP_X, J+1, I, NEIGHBOUR_R )    ! NEIGHBOUR ROOT ELEMENT NUMBER
        
        CALL FIND_RANK(NEIGHBOUR_R, TARGET_RANK)     ! ENIGHBOUR'S RANK
    
        CALL MPI_ISEND(SEND_BUFF, SR_COUNT, MPI_DOUBLE_PRECISION, &
                        TARGET_RANK, NEIGHBOUR_R, MPI_COMM_WORLD, &
                        REQUEST, IERROR)
        
        NUM_SEND = NUM_SEND + 1
        
        ARRAY_OF_REQUESTS_SEND(NUM_SEND) = REQUEST
        
        
    ENDIF


END SUBROUTINE NON_BLOCKING_EXCAHNGE_INTERFACE_Y



!-----------------------------------------------------------------------
!> INSIDE ONE ELEMENT USE NON-BLOCKING COMMUNICATION TO UPDATE MPI INTERFACES
!! Y DIRECTION
!-----------------------------------------------------------------------
SUBROUTINE NON_BLOCKING_EXCAHNGE_NFLUX_Y(FLAG1, FLAG2, LELEM_K, &
                                            NUM_SEND, &
                                            ARRAY_OF_REQUESTS_SEND, &
                                            SEND_BUFF, RECV_BUFF)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: LELEM_K  !< LOCAL ELEMENT NUMBER (START WITH 0)
    INTEGER :: RELEM_K  ! ROOT ELEMENT NUMBER (START WITH 0)
    INTEGER :: I, J     ! ELEMENT COORDINATES
    INTEGER :: NEIGHBOUR_R  ! NEIGHBOUR ELEMENT NUMBER (ROOT)
    INTEGER :: TARGET_RANK  ! NEIGHBOUR'S RANK
    INTEGER :: SR_COUNT   ! Number of elements in send/recv buffer 
    INTEGER :: REQUEST  ! Communication request (handle).  
    INTEGER :: IERROR
    INTEGER :: NUM_SEND ! COUND ISEND TIMES
    INTEGER :: ARRAY_OF_REQUESTS_SEND(LOCAL_ELEM_NUM) ! RECORD COMMINICATION REQUEST HANDLE
    INTEGER :: STATUS(MPI_STATUS_SIZE)

    LOGICAL, INTENT(IN) :: FLAG1 !< MPI BOUNDARY FLAG (FACE1/3)
    LOGICAL, INTENT(IN) :: FLAG2 !< MPI BOUNDARY FLAG (FACE2/4)
    
    DOUBLE PRECISION :: SEND_BUFF(0:NMAX, NUM_OF_EQUATION)  !< SEND BUFFER
    DOUBLE PRECISION :: RECV_BUFF(0:NMAX, NUM_OF_EQUATION)  !< RECV BUFFER
    
    SR_COUNT = (1 + NMAX) * NUM_OF_EQUATION
    
    CALL INDEX_LOCAL_TO_GLOBAL(RANK, LELEM_K, RELEM_K)  ! CONVERT ELEMENT LOCAL INDEX TO GLOBAL INDEX

    IF(FLAG2) THEN  ! RECV INFORMATION REMOTELY
    
        ! USE MPI WILDCARDS: MPI_ANY_TAG, MPI_ANY_SOURCE
        CALL MPI_RECV(RECV_BUFF, SR_COUNT, MPI_DOUBLE_PRECISION, &
                        MPI_ANY_SOURCE, RELEM_K, &
                        MPI_COMM_WORLD, STATUS, IERROR)
        
    ENDIF
    
    IF(FLAG1) THEN
    
        CALL d2xy ( EXP_X, RELEM_K, J, I )  ! ELEMENT COORDINATE
        
        CALL xy2d ( EXP_X, J-1, I, NEIGHBOUR_R )    ! NEIGHBOUR ROOT ELEMENT NUMBER
        
        CALL FIND_RANK(NEIGHBOUR_R, TARGET_RANK)     ! ENIGHBOUR'S RANK
    
        CALL MPI_ISEND(SEND_BUFF, SR_COUNT, MPI_DOUBLE_PRECISION, &
                        TARGET_RANK, NEIGHBOUR_R, MPI_COMM_WORLD, &
                        REQUEST, IERROR)
        
        NUM_SEND = NUM_SEND + 1
        
        ARRAY_OF_REQUESTS_SEND(NUM_SEND) = REQUEST
        
    ENDIF


END SUBROUTINE NON_BLOCKING_EXCAHNGE_NFLUX_Y

END MODULE MESSAGE_EXCHAGE
