!-----------------------------------------------------------------------
!> @brief
!> Adjustable parameters in the program
!-----------------------------------------------------------------------
MODULE PARAM

   USE MPI

    IMPLICIT NONE
    
    ! variables you must change =======================================
    
    !< INITIAL MESH FILE------------------------------------------------
    CHARACTER(LEN=*), PARAMETER :: MESHFILE = "16_elements.msh"    
    !-------------------------------------------------------------------
    
    ! SET POLYNOMIAL ORDER----------------------------------------------
    INTEGER :: N = 4   !< POLYNOMIAL DEGREE IN X DIRECTION
    INTEGER :: M = 4    !< POLYNOMIAL DEGREE IN Y DIRECTION
    INTEGER :: NMAX = 4   !< MAXIMUM POLYNOMIAL DEGREE IN X DIRECTION
    INTEGER :: MMAX = 4   !< MAXIMUM POLYNOMIAL DEGREE IN Y DIRECTION
    !-------------------------------------------------------------------
    
    ! DOMIAN BOUNDARY---------------------------------------------------
    DOUBLE PRECISION :: GX_L = 0.0D0     !< LEFT DOMAIN BOUNDARY
    DOUBLE PRECISION :: GX_R = 1.0D0     !< LEFT DOMAIN BOUNDARY
    
    DOUBLE PRECISION :: GY_L = 0.0D0     !< LEFT DOMAIN BOUNDARY
    DOUBLE PRECISION :: GY_R = 1.0D0     !< LEFT DOMAIN BOUNDARY
    !-------------------------------------------------------------------
    
    ! NUMBER OF ELEMENT-------------------------------------------------
    INTEGER :: EXP_X = 2       !< WRITE ELEMENT NUMBER IN X DIRECTION EXPONENTIAL ORDER, I.E., 2^(EXP_X)
    INTEGER :: EXP_Y = 2       !< WRITE ELEMENT NUMBER IN Y DIRECTION EXPONENTIAL ORDER, I.E., 2^(EXP_Y)
    !-------------------------------------------------------------------
    
    ! TIME--------------------------------------------------------------
!    DOUBLE PRECISION :: T_TOTAL = 0.0d0     !< TOTAL TIME INTEGRAL
!    INTEGER :: NT = 0                   !< TIME STEP NUMBER
!    DOUBLE PRECISION :: T_TOTAL = (2.0e-4)    !< TOTAL TIME INTEGRAL
!    INTEGER :: NT = 1                    !< TIME STEP NUMBER
    DOUBLE PRECISION :: T_TOTAL = 1.0D0     !< TOTAL TIME INTEGRAL
    INTEGER :: NT = 100000                   !< TIME STEP NUMBER
    !-------------------------------------------------------------------
    
    ! OUTPUT DATA FREQUENCY---------------------------------------------
    INTEGER :: OUTPUT_FREQUENCY = 1000    !< OUTPUT DATA FREQUENCY, SHOULD <= NT
    !-------------------------------------------------------------------
    
    !===================================================================
    
    ! variables you could change =======================================
    
    !< WAVE SPEED ------------------------------------------------------
    DOUBLE PRECISION, PARAMETER :: C=1.0D0  !< SOUND SPEED
    !-------------------------------------------------------------------
    
    ! SET ADAPATION ----------------------------------------------------
    INTEGER, PARAMETER :: SPLIT_MAX_NUM = 3     !< MAXIMUM SPLIT NUMBER
    !-------------------------------------------------------------------
    
    ! VERIFICATION -----------------------------------------------------
    LOGICAL :: VERIFICATION_SWITCH = .TRUE. !< VERIFICATE YOUR RESULTS. 
                                             !< NOTE: THIS SHOULD BE TURNED OFF IF YOU ARE USING REFLECTIVE BOUNDARY CONDITION
    !-------------------------------------------------------------------
    
    !===================================================================
    
    ! variables you do not need to adjust ==============================
    
    ! ------------------------------------------------------------------
    INTEGER :: NUM_OF_EQUATION = 3  !< NUMBER OF EQUATION
    !-------------------------------------------------------------------
    
    ! BOUNDARY FLAGS----------------------------------------------------
    INTEGER, PARAMETER :: NONE_L = 0
    INTEGER, PARAMETER :: NONE_R = 0
    !-------------------------------------------------------------------
    
    ! .MSH FILE PATH----------------------------------------------------
    CHARACTER(*), PARAMETER :: FILEPLACE = "./gmsh_files/"
    !-------------------------------------------------------------------
    
    ! TECPLOT FILE PATH-------------------------------------------------
    CHARACTER(*), PARAMETER :: OUTPUT_PLACE = "./output/"
    !-------------------------------------------------------------------
    
    ! SET MPI ----------------------------------------------------------
    INTEGER :: RANK     !< PROCESSOR RANK
    INTEGER :: NUM_PROC !< NUMBER OF PROCESSOR
    INTEGER :: IERROR 
    !-------------------------------------------------------------------
    
    !< OTHER PARAMETER --------------------------------------------------
    DOUBLE PRECISION, PARAMETER :: PI=4.0D0*DATAN(1.0D0)
    !-------------------------------------------------------------------
    
    !===================================================================

END MODULE PARAM
