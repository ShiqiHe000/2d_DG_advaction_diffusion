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

IMPLICIT NONE 

INTEGER :: FILE_NUM2 = 1 ! FILE NUMBER

CONTAINS

!-----------------------------------------------------------------------
!> INTERPOLATE THE SOLUTION TO 4 CORNERS OF THE CELL AND OUTPUT TO FILES
!-----------------------------------------------------------------------
SUBROUTINE WRITE_MESH_4(NEL_TOTAL, X_GLOBAL, Y_GLOBAL, PLEVELX, PLEVELY, &
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
    
    DOUBLE PRECISION :: SOLU_INT_L(0:NMAX, NUM_OF_EQUATION)
    DOUBLE PRECISION :: SOLU_INT_R(0:NMAX, NUM_OF_EQUATION)
    
    DOUBLE PRECISION :: X(4, NEL_TOTAL)  ! X COORDINATES OF FOUR NODES OF ONE ELEMENT  
    DOUBLE PRECISION :: Y(4, NEL_TOTAL)  ! Y COORDINATES OF FOUR NODES OF ONE ELEMENT  
    
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
    !-------------------------------------------------------------------

    SOLU_INT_L = 0.0D0; SOLU_INT_R = 0.0D0
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
    DO IEL = 0, NEL_TOTAL
        X(1, IEL+1) = X_GLOBAL(1, IEL)
        X(2, IEL+1) = X_GLOBAL(1, IEL)
        X(3, IEL+1) = X_GLOBAL(2, IEL)
        X(4, IEL+1) = X_GLOBAL(2, IEL)
    ENDDO
    
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
    
    !-------------------------------------------------------------------
    
    
    ! close file
    CALL h5fclose_f(file_id, error) ! file
    
    ! Close FORTRAN predefined datatypes.
     CALL h5close_f(error)


END SUBROUTINE WRITE_MESH_4

END MODULE OUTPUT_HDF5
