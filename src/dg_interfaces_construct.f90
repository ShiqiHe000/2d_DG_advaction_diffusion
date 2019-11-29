!-----------------------------------------------------------------------
!> @brief 
!> Construct interfaces. 
!! Incluse interior interfaces between elements and interfaces on the 
!! boundaries
!! Interpolate the solutions to boundaries, and apply boundary conditions.
!-----------------------------------------------------------------------

MODULE INTERFACES_CONSTRUCT

USE MPI
USE BASIS

IMPLICIT NONE 

CONTAINS

!-----------------------------------------------------------------------
!> Use Lagrange Polynomial Interplants to obtain the solution on the 
!! boundaries based on the interior collocation nodes.
!! Also enforce the boundary conditions on the boundary elements
!! This subroutine only construct interfaces in X direction.
!-----------------------------------------------------------------------
SUBROUTINE CONSTRUCT_INTERFACES_X(N1, M1, N_EQU, SOLUTION, &
                                LAG_1, LAG_2, SOLU_INT_L, SOLU_INT_R )

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N1    !< POLY ORDER OF CURRENT DIRECTION
    INTEGER, INTENT(IN) :: M1    !< POLY ORDER OF VERTICAL DIRECTION
    INTEGER, INTENT(IN) :: N_EQU   !< NUMBER OF EQUATIONS
    
    DOUBLE PRECISION :: SOLUTION(0:N1, 0:M1, N_EQU)   ! INTERIOR SOLUTIONS
    
    INTEGER :: S, J
    
    DOUBLE PRECISION :: LAG_1(0:N1) !< LAGRANGE INTERPOLATES TO THE LEFT BOUNDARY
    DOUBLE PRECISION :: LAG_2(0:N1) !< LAGRANGE INTERPOLATES TO THE RIGHT BOUNDARY
    
    DOUBLE PRECISION :: SOLU_INT_L(0:M1, N_EQU) ! SOLUTION ON THE LEFT BOUNDARY
    DOUBLE PRECISION :: SOLU_INT_R(0:M1, N_EQU) ! SOLUTION ON THE RIGHT BOUNDARY
    
    
    ! INTERPOLATE INTERIOR SOLUTION TO INTERFACES-----------------------
    DO S=1, N_EQU
    
        DO J=0, M1
            ! conventional vector vector multiplication
!            CALL INTERPOLATE_TO_BOUNDARY(N1, SOLUTION(:, J, S), &
!                                    LAG_1, SOLU_INT_L(J, S))
!            CALL INTERPOLATE_TO_BOUNDARY(N1, SOLUTION(:, J, S), &
!                                    LAG_2, SOLU_INT_R(J, S))
             
            ! use blas for vector vector multiplication 
            CALL INTERPOLATE_TO_BOUNDARY_BLAS(N1, SOLUTION(:, J, S), &
                                    LAG_1, SOLU_INT_L(J, S))
            CALL INTERPOLATE_TO_BOUNDARY_BLAS(N1, SOLUTION(:, J, S), &
                                    LAG_2, SOLU_INT_R(J, S))
        ENDDO
    ENDDO
    !-------------------------------------------------------------------

END SUBROUTINE CONSTRUCT_INTERFACES_X

!-----------------------------------------------------------------------
!> Use Lagrange Polynomial Interplants to obtain the solution on the 
!! boundaries based on the interiol collocation nodes.
!! Also enforce the boundary conditions on the boundary elements
!! This subroutine only construct interfaces in Y direction.
!-----------------------------------------------------------------------
SUBROUTINE CONSTRUCT_INTERFACES_Y(N1, M1, N_EQU, SOLUTION, &
                                LAG_1, LAG_2, SOLU_INT_L, SOLU_INT_R )

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N1    !< POLY ORDER OF CURRENT DIRECTION
    INTEGER, INTENT(IN) :: M1    !< POLY ORDER OF VERTICAL DIRECTION
    INTEGER, INTENT(IN) :: N_EQU   !< NUMBER OF EQUATIONS
    
    DOUBLE PRECISION :: SOLUTION(0:N1, 0:M1, N_EQU)   ! INTERIOR SOLUTIONS
    
    INTEGER :: S, I
    
    DOUBLE PRECISION :: LAG_1(0:M1) !< LAGRANGE INTERPOLATES TO THE LEFT BOUNDARY
    DOUBLE PRECISION :: LAG_2(0:M1) !< LAGRANGE INTERPOLATES TO THE RIGHT BOUNDARY
    
    DOUBLE PRECISION :: SOLU_INT_L(0:N1, N_EQU) ! SOLUTION ON THE LEFT BOUNDARY
    DOUBLE PRECISION :: SOLU_INT_R(0:N1, N_EQU) ! SOLUTION ON THE RIGHT BOUNDARY
    
    ! INTERPOLATE INTERIOR SOLUTION TO INTERFACES-----------------------
    DO S=1, N_EQU
    
        DO I=0, N1
            ! conventional vector vector multiplication
            CALL INTERPOLATE_TO_BOUNDARY(M1, SOLUTION(I, :, S), &
                                    LAG_1, SOLU_INT_L(I, S))
            CALL INTERPOLATE_TO_BOUNDARY(M1, SOLUTION(I, :, S), &
                                    LAG_2, SOLU_INT_R(I, S))
            
            ! use blas for vector vector multiplication 
            CALL INTERPOLATE_TO_BOUNDARY_BLAS(M1, SOLUTION(I, :, S), &
                                    LAG_1, SOLU_INT_L(I, S))
            CALL INTERPOLATE_TO_BOUNDARY_BLAS(M1, SOLUTION(I, :, S), &
                                    LAG_2, SOLU_INT_R(I, S))

        ENDDO
    ENDDO
    !-------------------------------------------------------------------

END SUBROUTINE CONSTRUCT_INTERFACES_Y



END MODULE INTERFACES_CONSTRUCT
