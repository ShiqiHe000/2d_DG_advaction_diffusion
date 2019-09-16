!-----------------------------------------------------------------------
!> @brief
!> Adjustable parameters in the program
!-----------------------------------------------------------------------
MODULE PARAM

    USE MPI

    IMPLICIT NONE
    
    !< INITIAL MESH FILE------------------------------------------------
    CHARACTER(LEN=*), PARAMETER :: MESHFILE = "two_boundary.msh"    
    !-------------------------------------------------------------------
    
    ! SET POLYNOMIAL ORDER----------------------------------------------
    INTEGER :: N = 14    !< POLYNOMIAL DEGREE IN X DIRECTION
    INTEGER :: M = 10    !< POLYNOMIAL DEGREE IN Y DIRECTION
    INTEGER :: MNMAX = 12   !< MAXIMUM POLYNOMIAL DEGREE
    !-------------------------------------------------------------------
    
    !< NUMBER OF EQUATION-----------------------------------------------
    INTEGER :: NUM_OF_EQUATION = 3
    !-------------------------------------------------------------------
    
    !< WAVE SPEED ------------------------------------------------------
    DOUBLE PRECISION, PARAMETER :: C=1.0D0
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
!    DOUBLE PRECISION :: K_X = DSQRT(2.0D0)/2.0D0    !< WAVE VECTOR IN X DIRECTION
!    DOUBLE PRECISION :: K_Y = DSQRT(2.0D0)/2.0D0    !< WAVE VECTOR IN Y DIRECTION
    DOUBLE PRECISION :: K_X = 1.0D0    !< WAVE VECTOR IN X DIRECTION
    DOUBLE PRECISION :: K_Y = 0.0D0    !< WAVE VECTOR IN Y DIRECTION
    !-------------------------------------------------------------------
    
    ! TIME--------------------------------------------------------------
!    DOUBLE PRECISION :: T_TOTAL = (2.6e-3)     !< TOTAL TIME INTEGRAL
!    INTEGER :: NT = 1                    !< TIME STEP NUMBER
    DOUBLE PRECISION :: T_TOTAL = 2.0D0     !< TOTAL TIME INTEGRAL
    INTEGER :: NT = 1000                   !< TIME STEP NUMBER
    !-------------------------------------------------------------------
    
    ! SET ADAPATION ----------------------------------------------------
    INTEGER, PARAMETER :: SPLIT_MAX_NUM = 3     !< MAXIMUM SPLIT NUMBER
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

END MODULE PARAM
