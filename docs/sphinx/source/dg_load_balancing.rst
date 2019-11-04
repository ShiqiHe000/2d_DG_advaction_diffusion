Dynamic Load-balancing 
***************************************

Motivation
======================
Load balance is one of the major challenges for the efficient supercomputer, especially for applications that exhibit workload variations.
Load imbalance is an considerable impedance on the path towards higher degree of parallelism. 

In particular, when load conditions changes dynamically, efficient mesh partitioning becomes an indispensible part of scalable design. 

Goals
=======================
Fluid dynamic application in the field of industrial engineering require high degrees of parallelism to achieve an acceptable time to solution for a large problem size. 
Typical mesh-based approaches therefore rely on suitable partitioning strategies to distribute the computational load across the set of processes. 

Therefore, an efficient load-balancing approach aims to achieve two goals:

- The work load should be distribute evenly 
        * avoid waiting times of processing units

- At the same time the interfacing boundaries between partitions should be as small as possible.
        * minimize the time spend in communication

The optimization problem is **NP-hard**.


Two popular approaches
================================================
Graph-based Algorithm
-------------------------------------------
A popular choice for graph-based partition is ParMetis_.

**ParMetis** performing very well for mesh partition for a long time. However, since ParMetis require global knowledge of the mesh, with an increassing number of processes, graph-based partitioning algorithms seem to reach their scalability limits. 
The memory consumption grows linearly with teh graph size, raising the need for alternatives which could avoid this problem. Such methods are based on `space-filling curves`_ (SFCs).

.. _ParMetis : http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview

.. _`space-filling curves` : https://en.wikipedia.org/wiki/Space-filling_curve


Space-filling curves (SFCs) based algorithm
----------------------------------------------
SFCs reduce the partitioning problem from n dimension to one dimension. 
The remaining tast, the so-called 1D patitioning problem or *chains-on-chains* partiitoning problem, is to decompose a 1D workload array into consecutive, balanced partitions. 


Advantages
^^^^^^^^^^^^^^^^^^^
- Good Locality
        * SFCs map the 1D unit interval onto a higher dimensional space such that neighboring points on the unit interval are also neighboring points in the target space. 

- Acceptable communication overhead
        * SFCs ignores the edges of full graph imformation. 
          It relies on the spatial properties of the curve to ensure a reasonable partition shape. Tirthapura et al. demonstrated that the upper limit of expected remote accesses in SFC partiiotned domains are acceptable :cite:`1`.  

- Low memory using
        * Taking the good locality of SFCs, the global imformation (full graph information) needed by Graph-based algorithm can be abandoned. Thus, SFCs opens a path towards low-memory partitioning strategies. 

Implementation to DG Solver
========================================
The numerical approximation of wave equation is a hp-adaptive approach. 
That is, elements can split or merge (h-adaptive) according to therequired resolution. Also, they can raise or decrease the polynomial degree (p-adaptive) to adjust the convergence rate. 

Due to the hp-adaptivity, different element can have differernt individual computation times, load imbalance is introuduce to this application. Re-meshing and domain partitioning are not avoidable. 

With the help of a SFC, the 2D structured mesh partitioning problem can be reduced to a 1D chains-on-chains partitioning (CCP) problem. Knowing the index od an element, its neighbours indices can be computed locally. 

Hilbert Curve
------------------------------
There are many SFCs, for example `Morton Curve`_ (z-curve) and `Hilbert Curve`_.

We choose Hilbert Curve as our SFC. Although Hilbert ordering is less efficient (with flip and rotation) than Morton Curve, Hilbert Curve brings out a better locality (no sudden "jump"). 

.. image:: /image/motorn_curve.png
.. image:: /image/Hilbert_curve.png

.. _`Morton Curve` : https://en.wikipedia.org/wiki/Z-order_curve

.. _`Hilbert Curve` : https://en.wikipedia.org/wiki/Hilbert_curve


References
===================================
.. bibliography:: dg_refs.bib
