!-----------------------------------------------------------------------
!> @brief
!> Adjustable parameters in the program
!-----------------------------------------------------------------------
MODULE PARAM

   USE MPI

    IMPLICIT NONE
    
    ! variables you could change =======================================
    
    !< INITIAL MESH FILE------------------------------------------------
    CHARACTER(LEN=*), PARAMETER :: MESHFILE = "two_boundary.msh"    
    !-------------------------------------------------------------------
    
    ! SET POLYNOMIAL ORDER----------------------------------------------
    INTEGER :: N = 8    !< POLYNOMIAL DEGREE IN X DIRECTION
    INTEGER :: M = 8    !< POLYNOMIAL DEGREE IN Y DIRECTION
    INTEGER :: MNMAX = 12   !< MAXIMUM POLYNOMIAL DEGREE
    !-------------------------------------------------------------------
    
    ! DOMIAN BOUNDARY---------------------------------------------------
    DOUBLE PRECISION :: GX_L = 1.0D0     !< LEFT DOMAIN BOUNDARY
    DOUBLE PRECISION :: GX_R = 2.0D0     !< LEFT DOMAIN BOUNDARY
    
    DOUBLE PRECISION :: GY_L = 1.0D0     !< LEFT DOMAIN BOUNDARY
    DOUBLE PRECISION :: GY_R = 2.0D0     !< LEFT DOMAIN BOUNDARY
    
!    DOUBLE PRECISION :: GX_L = -1.0D0     !< LEFT DOMAIN BOUNDARY
!    DOUBLE PRECISION :: GX_R =  1.0D0     !< LEFT DOMAIN BOUNDARY
    !-------------------------------------------------------------------
    
    ! NUMBER OF ELEMENT-------------------------------------------------
    INTEGER :: EXP_X = 1       !< WRITE ELEMENT NUMBER IN X DIRECTION EXPONENTIAL ORDER, I.E., 2^(EXP_X)
    INTEGER :: EXP_Y = 1       !< WRITE ELEMENT NUMBER IN Y DIRECTION EXPONENTIAL ORDER, I.E., 2^(EXP_Y)
    !-------------------------------------------------------------------
    
    !< WAVE SPEED ------------------------------------------------------
    DOUBLE PRECISION, PARAMETER :: C=1.0D0
    !-------------------------------------------------------------------
    
    ! TIME--------------------------------------------------------------
!    DOUBLE PRECISION :: T_TOTAL = 0.0d0     !< TOTAL TIME INTEGRAL
!    DOUBLE PRECISION :: T_TOTAL = (2.0e-4)     !< TOTAL TIME INTEGRAL
!    INTEGER :: NT = 1                    !< TIME STEP NUMBER
    DOUBLE PRECISION :: T_TOTAL = 2.0D0     !< TOTAL TIME INTEGRAL
    INTEGER :: NT = 10000                   !< TIME STEP NUMBER
    !-------------------------------------------------------------------
    
    ! SET ADAPATION ----------------------------------------------------
    INTEGER, PARAMETER :: SPLIT_MAX_NUM = 3     !< MAXIMUM SPLIT NUMBER
    !-------------------------------------------------------------------
    
    !===================================================================
    
    ! variables you do not need to adjust ==============================
    
    !< NUMBER OF EQUATION-----------------------------------------------
    INTEGER :: NUM_OF_EQUATION = 3
    !-------------------------------------------------------------------
    
    ! BOUNDARY FLAGS----------------------------------------------------
    INTEGER, PARAMETER :: NONE_L = 0
    INTEGER, PARAMETER :: NONE_R = 0
    !-------------------------------------------------------------------
    
    ! SET MPI ----------------------------------------------------------
    INTEGER :: RANK     !< PROCESSOR RANK
    INTEGER :: NUM_PROC !< NUMBER OF PROCESSOR
    INTEGER :: IERROR 
    !-------------------------------------------------------------------
    
    !< OUTPUT FILE ------------------------------------------------------
    INTEGER :: FRAME=1  
    !-------------------------------------------------------------------
    
    !< OTHER PARAMETER --------------------------------------------------
    DOUBLE PRECISION, PARAMETER :: PI=4.0D0*DATAN(1.0D0)
    !-------------------------------------------------------------------
    
    !===================================================================

END MODULE PARAM
