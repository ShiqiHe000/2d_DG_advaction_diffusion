!-----------------------------------------------------------------------
!> @brief
!> Use HDF5 format to output data.
!! Store each variables in a dataset
!-----------------------------------------------------------------------

MODULE OUTPUT_HDF5

USE HDF5
USE MPI
USE PARAM
USE LOCAL_STORAGE
USE POLY_LEVEL_AND_ORDER
USE NODAL_2D_STORAGE

IMPLICIT NONE 

INTEGER :: FILE_NUM2 = 1 ! FILE NUMBER

! HDF5--------------------------------------------------------------
INTEGER :: ERROR
INTEGER(HID_T) :: file_id ! File identifier 
INTEGER(HID_T) :: prp_id  ! File Property list identifier 
INTEGER(HID_T) :: dataspace_id  ! Dataspace identifier 
INTEGER(HID_T) :: dset_id ! Dataset identifier
INTEGER(HID_T) :: memspace ! memory identifier 
    
INTEGER(HSIZE_T), DIMENSION(2) :: data_dim  ! TOTAL DATA DIMENSSION FOR THE DATASET
INTEGER(HSIZE_T), DIMENSION(2) :: OFFSET    ! position to write
INTEGER(HSIZE_T), DIMENSION(2) :: local_data_dim    ! local data dimension

CHARACTER(LEN=6), PARAMETER :: X_DATASET = "xcoord"
CHARACTER(LEN=6), PARAMETER :: Y_DATASET = "ycoord"
CHARACTER(LEN=8), PARAMETER :: P_DATASET = "pressure"
CHARACTER(LEN=1), PARAMETER :: U_DATASET = "u"
CHARACTER(LEN=1), PARAMETER :: V_DATASET = "v"
!-------------------------------------------------------------------


CONTAINS

!-----------------------------------------------------------------------
!> INTERPOLATE THE SOLUTION TO 4 CORNERS OF THE CELL AND OUTPUT TO FILES
!-----------------------------------------------------------------------
SUBROUTINE WRITE_MESH_P4(NEL_TOTAL, X_GLOBAL, Y_GLOBAL, PLEVELX, PLEVELY, &
                        SOLUTION_ALL, T)
    
    IMPLICIT NONE 
    
    CHARACTER(LEN=15) :: FILENAME
    
    INTEGER :: ELEM, IEL
    INTEGER :: PORDERX, PORDERY
    INTEGER, INTENT(IN) :: NEL_TOTAL    !< TOTAL NUMBER OF ELEMENT
    
    INTEGER :: PLEVELX(0:NEL_TOTAL-1)  !< POLYNOMIAL LEVEL IN X(CAN BE TRANSFORM TO POLY ORDER)
    INTEGER :: PLEVELY(0:NEL_TOTAL-1)  !< POLYNOMIAL LEVEL IN Y(CAN BE TRANSFORM TO POLY ORDER)
    
    
    DOUBLE PRECISION :: X_GLOBAL(2, 0:NEL_TOTAL-1)  !< ELEMENT X COORDINATES
    DOUBLE PRECISION :: Y_GLOBAL(2, 0:NEL_TOTAL-1)  !< ELEMENT Y COORDINATES
    
    DOUBLE PRECISION :: SOLUTION_ALL(0:NMAX, 0:MMAX, NUM_OF_EQUATION, 0:NEL_TOTAL-1)    !< NUMERICAL APPROXIMATION
    
    DOUBLE PRECISION :: T !< CURRENT TIME
    
    DOUBLE PRECISION :: X(4, NEL_TOTAL)  ! X COORDINATES OF FOUR NODES OF ONE ELEMENT  
    DOUBLE PRECISION :: Y(4, NEL_TOTAL)  ! Y COORDINATES OF FOUR NODES OF ONE ELEMENT  
    
    DOUBLE PRECISION :: P(4, NEL_TOTAL) ! PRESSURE
    DOUBLE PRECISION :: U(4, NEL_TOTAL) ! X VELOCITY COMPONENT 
    DOUBLE PRECISION :: V(4, NEL_TOTAL) ! Y VELOCITY COMPONENT
    
    
!    SOLU_INT_L = 0.0D0; SOLU_INT_R = 0.0D0
    X = 0.0D0; Y = 0.0D0
    
    ! create file name--------------------------------------------------
    WRITE(FILENAME,FMT='(''aoutput'',I5.5,''.h5'')') FILE_NUM2
    !-------------------------------------------------------------------
    
    ! Initialize FORTRAN interface--------------------------------------
    CALL h5open_f(IERROR)
    !-------------------------------------------------------------------
    
    ! create file property list-----------------------------------------
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, prp_id, error)  
    CALL h5pset_fapl_mpio_f(prp_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)   ! let MPI have access
    !-------------------------------------------------------------------
    
    ! create file ------------------------------------------------------
    CALL h5fcreate_f(FILENAME, H5F_ACC_TRUNC_F, file_id, error, &  
                       access_prp = prp_id)
    CALL h5pclose_f(prp_id, error)
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    data_dim(1) = 4; data_dim(2) = ELEM_RANGE(NUM_PROC) + 1
    local_data_dim(1) = 4; local_data_dim(2) = NEL_TOTAL 
    !-------------------------------------------------------------------
    
    ! OUTPUT X COORDINATES----------------------------------------------
        
    ! PREPARE THE X COORD TO WRITE
    DO IEL = 0, NEL_TOTAL - 1
        X(1, IEL+1) = X_GLOBAL(1, IEL)
        X(2, IEL+1) = X_GLOBAL(1, IEL)
        X(3, IEL+1) = X_GLOBAL(2, IEL)
        X(4, IEL+1) = X_GLOBAL(2, IEL)
    ENDDO
    
    CALL WRITE_X_COORDS(NEL_TOTAL, X)
    !-------------------------------------------------------------------
    
    ! OUTPUT Y COORDINATES----------------------------------------------
        
    ! PREPARE THE X COORD TO WRITE
    DO IEL = 0, NEL_TOTAL - 1
        Y(1, IEL+1) = Y_GLOBAL(1, IEL)
        Y(2, IEL+1) = Y_GLOBAL(2, IEL)
        Y(3, IEL+1) = Y_GLOBAL(2, IEL)
        Y(4, IEL+1) = Y_GLOBAL(1, IEL)
    ENDDO
    
    CALL WRITE_Y_COORDS(NEL_TOTAL, Y)
    !-------------------------------------------------------------------
    
    ! OUTPUT SOLUTIONS -------------------------------------------------
    CALL CONSTRUCT_CORNERS(NEL_TOTAL, SOLUTION_ALL, PLEVELX, PLEVELY, P, U, V)
    
    CALL WRITE_PUV(NEL_TOTAL, P, P_DATASET)
    
    CALL WRITE_PUV(NEL_TOTAL, U, U_DATASET)
    
    CALL WRITE_PUV(NEL_TOTAL, V, V_DATASET)
    !-------------------------------------------------------------------
    
    ! close file
    CALL h5fclose_f(file_id, error) ! file
    
    ! Close FORTRAN predefined datatypes.
    CALL h5close_f(error)
    
    FILE_NUM2 = FILE_NUM2 + 1


END SUBROUTINE WRITE_MESH_P4

!-----------------------------------------------------------------------
!> WRITE X COORDINATES.
!! FOUR COLUMNS: 4 CORNER NODES
!! EACH ROW: ONE ELEMENT
!-----------------------------------------------------------------------
SUBROUTINE WRITE_X_COORDS(NEL_TOTAL, X)
    
    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: NEL_TOTAL
    
    DOUBLE PRECISION :: X(4, NEL_TOTAL)  ! X COORDINATES OF FOUR NODES OF ONE ELEMENT  

    ! create dataspace
    CALL h5screate_simple_f(2, data_dim, dataspace_id, error)
    
    ! create dataset (Creates a new dataset and links it to a location in the file. )
    CALL h5dcreate_f(file_id, X_DATASET, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, & 
                    error)
    CALL h5sclose_f(dataspace_id, error) 
    
    ! defines the location of write 
    OFFSET(1) = 0
    OFFSET(2) = (ELEM_RANGE(RANK) + 1) 
    
    ! Select hyperslab in the file.
    CALL h5dget_space_f(dset_id, dataspace_id, error)
    CALL h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, OFFSET, local_data_dim, error) 
    
    ! memory dataspace
    CALL h5screate_simple_f(2, local_data_dim, memspace, error)
    
    ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, prp_id, error) 
    CALL h5pset_dxpl_mpio_f(prp_id, H5FD_MPIO_COLLECTIVE_F, error)   ! Sets data transfer mode. 
    
    
    ! Writes raw data from a buffer to a dataset. 
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, X, local_data_dim, error, &  ! the dataset_dim here is the local dim (global dim also works)
                      & mem_space_id = memspace, file_space_id = dataspace_id, xfer_prp = prp_id)

    ! close resource
    CALL h5sclose_f(dataspace_id, error)   ! space
    CALL h5sclose_f(memspace, error)
    
    CALL h5pclose_f(prp_id, error)  ! property
    
    CALL h5dclose_f(dset_id, error)   ! dataset



END SUBROUTINE WRITE_X_COORDS


!-----------------------------------------------------------------------
!> WRITE Y COORDINATES.
!! FOUR COLUMNS: 4 CORNER NODES
!! EACH ROW: ONE ELEMENT
!-----------------------------------------------------------------------
SUBROUTINE WRITE_Y_COORDS(NEL_TOTAL, Y)
    
    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: NEL_TOTAL
    
    DOUBLE PRECISION :: Y(4, NEL_TOTAL)  ! X COORDINATES OF FOUR NODES OF ONE ELEMENT  

    ! create dataspace
    CALL h5screate_simple_f(2, data_dim, dataspace_id, error)
    
    ! create dataset (Creates a new dataset and links it to a location in the file. )
    CALL h5dcreate_f(file_id, Y_DATASET, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, & 
                    error)
    CALL h5sclose_f(dataspace_id, error) 
    
    ! defines the location of write 
    OFFSET(1) = 0
    OFFSET(2) = (ELEM_RANGE(RANK) + 1) 
    
    ! Select hyperslab in the file.
    CALL h5dget_space_f(dset_id, dataspace_id, error)
    CALL h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, OFFSET, local_data_dim, error) 
    
    ! memory dataspace
    CALL h5screate_simple_f(2, local_data_dim, memspace, error)
    
    ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, prp_id, error) 
    CALL h5pset_dxpl_mpio_f(prp_id, H5FD_MPIO_COLLECTIVE_F, error)   ! Sets data transfer mode. 
    
    
    ! Writes raw data from a buffer to a dataset. 
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Y, local_data_dim, error, &  ! the dataset_dim here is the local dim (global dim also works)
                      & mem_space_id = memspace, file_space_id = dataspace_id, xfer_prp = prp_id)

    ! close resource
    CALL h5sclose_f(dataspace_id, error)   ! space
    CALL h5sclose_f(memspace, error)
    
    CALL h5pclose_f(prp_id, error)  ! property
    
    CALL h5dclose_f(dset_id, error)   ! dataset



END SUBROUTINE WRITE_Y_COORDS



!-----------------------------------------------------------------------
!> INTERFACE ILLUSTRATION:
!!           X_R
!!      -------------
!!      |            |
!!      |            |
!! Y_L  |            | Y_R
!!      |            |
!!      |            |
!!      -------------
!!           X_L   
!-----------------------------------------------------------------------

SUBROUTINE CONSTRUCT_CORNERS(NEL_LOCAL, Q, PLEVELX, PLEVELY, P, U, V)

    USE INTERFACES_CONSTRUCT

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NEL_LOCAL    !< LOCAL ELEMENT NUMBER
    INTEGER :: PORDERX, PORDERY
    INTEGER :: IEL, K
      
    INTEGER, INTENT(IN) :: PLEVELX(0:NEL_LOCAL-1)  !< POLYNOMIAL LEVEL IN X(CAN BE TRANSFORM TO POLY ORDER)
    INTEGER, INTENT(IN) :: PLEVELY(0:NEL_LOCAL-1)  !< POLYNOMIAL LEVEL IN Y(CAN BE TRANSFORM TO POLY ORDER)
    
    DOUBLE PRECISION, INTENT(IN) :: Q(0:NMAX, 0:MMAX, NUM_OF_EQUATION, 0:NEL_LOCAL-1)   ! SOLUTIONS
    
    DOUBLE PRECISION :: SOLUTION_X_L(0:MMAX, NUM_OF_EQUATION)
    DOUBLE PRECISION :: SOLUTION_X_R(0:MMAX, NUM_OF_EQUATION)
    
        
    DOUBLE PRECISION :: P(4, NEL_LOCAL) ! PRESSURE
    DOUBLE PRECISION :: U(4, NEL_LOCAL) ! X VELOCITY COMPONENT 
    DOUBLE PRECISION :: V(4, NEL_LOCAL) ! Y VELOCITY COMPONENT
    
    
    DO IEL = 0, NEL_LOCAL -1
        
        CALL POLY_LEVEL_TO_ORDER(N, PLEVELX(IEL), PORDERX)
        CALL POLY_LEVEL_TO_ORDER(M, PLEVELY(IEL), PORDERY)

        CALL CONSTRUCT_INTERFACES_X(PORDERX, PORDERY, &
                                    NUM_OF_EQUATION, &
                                    Q(0:PORDERX, 0:PORDERY, :, IEL), &
                                    LAGRANGE_LEFT_T(0:PORDERX, PLEVELX(IEL)),&
                                    LAGRANGE_RIGHT_T(0:PORDERX, PLEVELX(IEL)), &
                                    SOLUTION_X_L(0:PORDERY, :), &
                                    SOLUTION_X_R(0:PORDERY, :) )
        
        ! PREESURE------------------------------------------------------
        CALL FOUR_CORNERS(PORDERY, PLEVELY(IEL), SOLUTION_X_L(0:PORDERY, 1), &
                            SOLUTION_X_R(0:PORDERY, 1), P(:, IEL+1))
        !---------------------------------------------------------------
        
        ! VELOCITY U----------------------------------------------------
        CALL FOUR_CORNERS(PORDERY, PLEVELY(IEL), SOLUTION_X_L(0:PORDERY, 2), &
                            SOLUTION_X_R(0:PORDERY, 2), U(:, IEL+1))
        !---------------------------------------------------------------
        
        ! VELOCITY V----------------------------------------------------
        CALL FOUR_CORNERS(PORDERY, PLEVELY(IEL), SOLUTION_X_L(0:PORDERY, 3), &
                            SOLUTION_X_R(0:PORDERY, 3), V(:, IEL+1))
        !---------------------------------------------------------------
        
    ENDDO

END SUBROUTINE CONSTRUCT_CORNERS


SUBROUTINE FOUR_CORNERS(PORDERY, LEVELY, QL, QR, Q_OUT)
    
    USE BASIS
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: PORDERY  !< POLY ORDER IN Y DIRECTION
    INTEGER, INTENT(IN) :: LEVELY  !< ELEMENT LEVEL IN Y DIRECTION
    
    DOUBLE PRECISION, INTENT(IN) :: QL(0:PORDERY)    !< SOLUTION ON THE Y LEFT INTERFACES
    DOUBLE PRECISION, INTENT(IN) :: QR(0:PORDERY)    !< SOLUTION ON THE Y RIGHT INTERFACES
    
    DOUBLE PRECISION :: Q_OUT(4)    !< VALUES ON THE FOUR CORNERS
    
    ! CORNER 1 & 2
    CALL INTERPOLATE_TO_BOUNDARY(PORDERY, QL, &
                                LAGRANGE_DOWN_T(0:PORDERY, LEVELY), Q_OUT(1))
    CALL INTERPOLATE_TO_BOUNDARY(PORDERY, QL, &
                                LAGRANGE_UP_T(0:PORDERY, LEVELY), Q_OUT(2))
    
    ! CORNER 3 & 4
    CALL INTERPOLATE_TO_BOUNDARY(PORDERY, QR, &
                                LAGRANGE_DOWN_T(0:PORDERY, LEVELY), Q_OUT(3))
    CALL INTERPOLATE_TO_BOUNDARY(PORDERY, QR, &
                                LAGRANGE_UP_T(0:PORDERY, LEVELY), Q_OUT(3))
    

END SUBROUTINE FOUR_CORNERS

!-----------------------------------------------------------------------
!> WRITE P/U/V 
!! FOUR COLUMNS: 4 CORNER NODES
!! EACH ROW: ONE ELEMENT
!-----------------------------------------------------------------------
SUBROUTINE WRITE_PUV(NEL_TOTAL, H, H_DATASET)
    
    IMPLICIT NONE 
    
    INTEGER, INTENT(IN) :: NEL_TOTAL    !< LOCAL ELEMENT NUMEBR
    
    DOUBLE PRECISION :: H(4, NEL_TOTAL)  !< PRESSURE/VELOCITY 
    
    CHARACTER(LEN=*) :: H_DATASET

    ! create dataspace
    CALL h5screate_simple_f(2, data_dim, dataspace_id, error)
    
    ! create dataset (Creates a new dataset and links it to a location in the file. )
    CALL h5dcreate_f(file_id, H_DATASET, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, & 
                    error)
    CALL h5sclose_f(dataspace_id, error) 
    
    ! defines the location of write 
    OFFSET(1) = 0
    OFFSET(2) = (ELEM_RANGE(RANK) + 1) 
    
    ! Select hyperslab in the file.
    CALL h5dget_space_f(dset_id, dataspace_id, error)
    CALL h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, OFFSET, local_data_dim, error) 
    
    ! memory dataspace
    CALL h5screate_simple_f(2, local_data_dim, memspace, error)
    
    ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, prp_id, error) 
    CALL h5pset_dxpl_mpio_f(prp_id, H5FD_MPIO_COLLECTIVE_F, error)   ! Sets data transfer mode. 
    
    
    ! Writes raw data from a buffer to a dataset. 
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, H, local_data_dim, error, &  ! the dataset_dim here is the local dim (global dim also works)
                      & mem_space_id = memspace, file_space_id = dataspace_id, xfer_prp = prp_id)

    ! close resource
    CALL h5sclose_f(dataspace_id, error)   ! space
    CALL h5sclose_f(memspace, error)
    
    CALL h5pclose_f(prp_id, error)  ! property
    
    CALL h5dclose_f(dset_id, error)   ! dataset



END SUBROUTINE WRITE_PUV

END MODULE OUTPUT_HDF5
