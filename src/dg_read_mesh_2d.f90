!-----------------------------------------------------------------------
!> @brief
!> Read_mesh module read mesh information directly from GMSH .msh file.
!! NOTE: the .msh file has to be in version 2 ASII format.  
!-----------------------------------------------------------------------

MODULE READ_MESH

USE MPI

CONTAINS

SUBROUTINE READ_MESH_2D
    
    USE PARAM, ONLY: MESHFILE
    
    IMPLICIT NONE 
    
    

END SUBROUTINE READ_MESH_2D


END MODULE READ_MESH
