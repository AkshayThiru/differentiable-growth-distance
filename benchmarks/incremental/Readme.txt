The code for the incremental growth distance solver was obtained from
https://sites.google.com/view/ongcjhome/download.


Fast Growth Distance Solver
===========================

1.	Introduction
The fast Growth Distance solver is based on "An Incremental Version
of Growth Distance" by C.J. Ong and Eugene Huang, proceedings of IEEE
conference on Robotics and Automation 1998 with subsequent modification
by Y. Zhang.  This code is a descent based algorithm.  It finds the
Growth Function by decreasing sigma via a series of steps using the 
faces of the objects and their adjacency.


2.	Usage

Files contained in this package:
Demo.c            -   Demonstration code
Interface.c       -   Interface routine to Growth Distance solver
Interface.h       -   Header file for interface routine
binarySearch.c    -   Binary search algorithm used by Growth Distance
                      solver
binarySearch.h    -   Header file for binary search code
Gamma.c           -   Data management code used by Growth Distance solver
Gamma.h           -   Header file for "Gamma" code
incrCode4c.c      -   Code used by incremental Growth Distance solver
incrCode4c.h      -   Header file for incremental code
Solver4c.c        -   Main Growth Distance solver routine
Solver4c.h        -   Header file for main Growth Distance solver code
Tolerance.h       -   Tolerance used by the above routines
Readme.txt        -   This file


The routine the user should call is "gd" in "interface.c". However, "gd"
expects object data in a fixed format.  For users who are familiar with
this code, it is possible to call "gd" directly with the exact data format.
It is also possible to generate the correct data format using
the routine "setup_gamma" from "interface.c".  A description of other
required input parameters of "gd" is given in the preamble of function "gd"
in "interface.c".  For each object, "setup_gamma" is invoked only once.
Please read documentation on "setup_gamma" in "interface.c" for more details.
All other routines are used internally by the Growth Distance solver.
An example of usage of this code is given in "Demo.c". Users are advised
to go through this file. 

This current version does not include the hierarchical division of complex
faces. As long as objects have fewer than 50 neighboring faces for each face, 
the performance of this code does not degrade significantly.

We have tested this code quite extensively on a family of objects and
 believe that it is bug-free. If you do find problems with this code,
please send details to mpeongcj@nus.sg

3.	Compiling the demonstration code
For testing, there is a "Demo.c" test program that demonstrates the
use of the function "gd".  A makefile is provided to simplify the
compilation of the demostration program.
