!-----------------------------------------------------------------------
!> @brief
!> Processor's local storages.
!! Generate the basis information (GL PIOINTS, WEIGHTS, DERIVATIVE 
!! MATRIX, etc.) for all polynomial level at the beginning. 
!! Therefore, after the adaptation process starts, we do not need to 
!! re-generate the basis, or store duplicate data.
!-----------------------------------------------------------------------

MODULE BASIS_STORAGE

USE MPI
USE NODAL_2D_STORAGE
USE PARAM, ONLY: N, M, NMAX, MMAX
USE POLY_LEVEL_AND_ORDER

IMPLICIT NONE

CONTAINS

SUBROUTINE CONSTRUCT_BASIS_STORAGE

    IMPLICIT NONE
    
    INTEGER :: I, J, K
    INTEGER :: PORDER
    
    ! GENERATE MAX POLY ORDER
    LEVEL_MAX_X = (NMAX - N) / 2 + 1
    LEVEL_MAX_Y = (MMAX - M) / 2 + 1
    
    
    ! BASIS TABLE GL_POINT----------------------------------------------
    ALLOCATE(GL_POINT_X_T(0:NMAX, LEVEL_MAX_X))
    ALLOCATE(GL_POINT_Y_T(0:MMAX, LEVEL_MAX_Y))
    
    GL_POINT_X_T = 0.0D0; GL_POINT_Y_T = 0.0D0
    
    ALLOCATE(GL_W_X_T(0:NMAX, LEVEL_MAX_X))
    ALLOCATE(GL_W_Y_T(0:MMAX, LEVEL_MAX_Y))
    
    GL_W_X_T = 0.0D0; GL_W_Y_T = 0.0D0
    
    ALLOCATE(FIRST_DER_X_T(0:NMAX, 0:NMAX, LEVEL_MAX_X))
    ALLOCATE(FIRST_DER_Y_T(0:MMAX, 0:MMAX, LEVEL_MAX_Y))
    
    FIRST_DER_X_T = 0.0D0; FIRST_DER_Y_T = 0.0D0
    
    ALLOCATE(M_FIRST_DER_X_T(0:NMAX, 0:NMAX, LEVEL_MAX_X))
    ALLOCATE(M_FIRST_DER_Y_T(0:MMAX, 0:MMAX, LEVEL_MAX_Y))
    
    M_FIRST_DER_X_T = 0.0D0; M_FIRST_DER_Y_T = 0.0D0
    
    ALLOCATE(LAGRANGE_LEFT_T(0:NMAX, LEVEL_MAX_X))
    ALLOCATE(LAGRANGE_RIGHT_T(0:NMAX, LEVEL_MAX_X))
    
    LAGRANGE_LEFT_T = 0.0D0; LAGRANGE_RIGHT_T = 0.0D0
    
    ALLOCATE(LAGRANGE_DOWN_T(0:MMAX, LEVEL_MAX_Y))
    ALLOCATE(LAGRANGE_UP_T(0:MMAX, LEVEL_MAX_Y))
    
    LAGRANGE_DOWN_T = 0.0D0; LAGRANGE_UP_T = 0.0D0
    !-------------------------------------------------------------------
    
    ! X DIRECTION-------------------------------------------------------
    DO K=1, LEVEL_MAX_X
    
        CALL POLY_LEVEL_TO_ORDER(N, K, PORDER)
        
        CALL GET_NODAL_2D_STORAGE_BASIS(PORDER, &
                                        GL_POINT_X_T(0:PORDER, K), &
                                        GL_W_X_T(0:PORDER, K), &
                                        FIRST_DER_X_T(0:PORDER, 0:PORDER, K))
        
        CALL GET_NODAL_2D_STORAGE_EXTENDS(PORDER, &
                                          LAGRANGE_LEFT_T(0:PORDER, K), &
                                          LAGRANGE_RIGHT_T(0:PORDER, K), &
                                          GL_POINT_X_T(0:PORDER, K))
                                          
        DO J=0, PORDER
            DO I=0, PORDER
                M_FIRST_DER_X_T(I, J, K) = &
                - FIRST_DER_X_T(J, I, K) * GL_W_X_T(J, K) / GL_W_X_T(I, K)
            ENDDO
        ENDDO
        
    ENDDO
    !-------------------------------------------------------------------
    
    ! Y DIRECTION-------------------------------------------------------
    DO K=1, LEVEL_MAX_Y
    
        CALL POLY_LEVEL_TO_ORDER(M, K, PORDER)
        
        CALL GET_NODAL_2D_STORAGE_BASIS(PORDER, &
                                        GL_POINT_Y_T(0:PORDER, K), &
                                        GL_W_Y_T(0:PORDER, K), &
                                        FIRST_DER_Y_T(0:PORDER, 0:PORDER, K))
        
        CALL GET_NODAL_2D_STORAGE_EXTENDS(PORDER, &
                                          LAGRANGE_DOWN_T(0:PORDER, K), &
                                          LAGRANGE_UP_T(0:PORDER, K), &
                                          GL_POINT_Y_T(0:PORDER, K))
                                          
        DO J=0, PORDER
            DO I=0, PORDER
                M_FIRST_DER_Y_T(I, J, K) = &
                - FIRST_DER_Y_T(J, I, K) * GL_W_Y_T(J, K) / GL_W_Y_T(I, K)
            ENDDO
        ENDDO
    ENDDO
    !-------------------------------------------------------------------
    
    
    DEALLOCATE(FIRST_DER_X_T, FIRST_DER_Y_T)


END SUBROUTINE CONSTRUCT_BASIS_STORAGE


END MODULE BASIS_STORAGE
