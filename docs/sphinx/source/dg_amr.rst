AMR Strategies
**********************************

Introduction
===========================
Why AMR?
-----------------------
In order to effectively utilize the computational resources while remaining the flexibility in solving complex geometries and the prescribed accuracy, **Adaptive Mesh Refinement (AMR)** is invoked to focus the computational effort and memory usage to where it is needed.

Three main algorithms
---------------------------------
Three main algorithms have emerged overtime, which we can call them: **unstructured (U), block-structured (s)**, and hierarchical or **tree-based (T)** AMR.

UAMR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Unstructured mesh. Traditionally use graph-based partitioning algorithm, now are supplementing by fast algorithms based on coordinate partitioning and SFCs. 


SAMR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|pic1|  |pic2|

.. |pic1| image:: /image/AMR/block1.png
        :width: 45%

.. |pic2| image:: /image/AMR/block2.png
        :width: 45%

A sequence of nested structured grids at different hierachies or levels are overlapped with or patched onto each other. 

A tree-like data structure is used to facilitate the communication (transfer information) between the regular Cartesian grids at the various hierachies.
Each node in this tree-like data structure represents an entire grid rather than simply a cell. 

.. image:: /image/AMR/block_intergrate_patching.png

Pros:
        * Each node in the tree structure represents an entire grid enables the solver to solve the structured grids efficiently.  

Cons:
        * Communication patterns between levels can be complex.
        * Algorithm complexity can be substantial.
        * Due to the clustering methods used to define the sub-grids, 
          portions of the cumputational dmain covered by a highly refined mesh when it is not needed, resulting a wasted computational effort.

Library
+++++++++++++++
* Chombo
* PARAMESH
* SAMRAI

TAMR
^^^^^^^^^^^^^^^^^^^^^^^^^^

|pic3|  |pic4|

.. |pic3| image:: /image/AMR/tree1.png
        :width: 45%

.. |pic4| image:: /image/AMR/tree2.png
        :width: 45%

A qurad-tree/oct-tree data structure is used in 2D/3D to represent the grid hierarchies. Each node stands for a individual cell. 

.. image:: /image/AMR/quadtree_illustration.gif

Pros:
        * Mesh can be locally refined (increase storage savings)
        * Better control of the grid resolution (comparing with SAMR)

Cons:
        * In conventional quard-tree/oct-tree discretization, the connectivity information between individual cell and its neighbours needs to be stored explicitly. (oct-tree each cell 19 words of computer memory)

        * large memory overhead to maintain tree-data structures. 

        * Difficult to parallelize.
                - data moving: distruct and rebuild the linker.
                - neighbour finding: need to traverse the tree to locate the closet ancestor (what if ancestor is on another processor?).  

Library
+++++++++++++++++++++
* p4est
* Zoltan

Data structure: Quardtree/Octree
=======================================
Data structure Classifications
-----------------------------------
.. image:: /image/AMR/data_structure_Classification.jpg

Quardtree definition
--------------------------------
A `quadtree`_ is a tree data structure in which each internal node has exactly four children. Quadtrees are the two-dimensional analog of octrees and are most often used to partition a two-dimensional space by recursively subdividing it into four quadrants or regions.

.. _`quadtree` : https://en.wikipedia.org/wiki/Quadtree


Tree-based AMR algorithm 
=========================================
Objectives
--------------------------------------------
* Reduce the memory overhead required to maintain the information embodies in the tree structure.
* Rapid and easy access to the information stored in the tree. 

p4est
-----------------------------------
Linear octree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. image:: /image/AMR/p4est_linear_tree.png

Only store the leaves of the octree ("linear" octree).

.. image:: /image/AMR/schematic_for_refinement_p4est.png

Full Threaded Tree (FTT)
-----------------------------------
Memory requirement: :math:`2\frac{3}{8}` words per cell (conventional 19 words per cell).

The actual number of traversed levels required to find a neighbour never exceeds one. 

.. image:: /image/AMR/FTT_Oct.png
 
Cell-Based Structured Adaptive Mesh Refinement
-------------------------------------------------
Optimized FTT.

Cartesian-like indices are used to identify each cell. With these stored indices, the information on the parent, children and neighbours of a given cell can be accessed simply and efficiently.

Memory requirement: :math:`\frac{5}{8}` words per cell. 

.. image:: /image/AMR/CSAMR_Oct.png

Octant coordinate calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The indices of the four children octs :math:`(i_s, j_s)`

.. math::
        (i_s, j_s) = \left \{ (2i+m, 2j+n)| m = 0, 1; n = 0, 1 \right \}

The parent of a oct :math:`(i_p, j_p)`

.. math::
        (i_p, j_p) = \left ( int[\frac{i}{2}], int[\frac{j}{2}] \right )


Neighbour finding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. image:: /image/AMR/three_circumstances_of_neighbour.png

Cell3 find east neighbour:

(1). (i+1, j) -- hash table -- cell exsit (Y/N)?

(2). If Yes. 
        * Neighbour is the Northwest (NW) cell of cell(i+1, j) -- if this cell is a leaf (Y/N)?

                - Yes -- over

                - No -- two neighbours (NW, SW)


(3). If No. 
        * Neighbour cell has a larger size. cell number is :math:`\left ( int[\frac{i+1}{2}], int[\frac{j}{2}] \right )`


At most, two search are sufficient to find a neighbour of a give cell. Half of the neighbours can be reached without consulting the hash table. Statistically, the average number of searches required to find a neighbour of a given cell is one. 
