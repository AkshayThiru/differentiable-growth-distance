#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Interface.h"

/* This is a sample program for demonstrating the kind of data and
   initiliazations required by the Growth Distance 3D Distance
   Computation Routine. */

/* Value of PI/180. */
#define PIOVER180   0.017453292519943



/* =============================================================== */
/* In this demonstration, convex objects are declared as intersection
   of half-spaces.  The specification of a half-space is given
   by the outward normal vector and a corresponding point on
   the plane.

   Data declarations: All variables needed to successfully run the
   Growth Distance solver are defined here.

   DATA REQUIRED
   =============

   faces1,
   faces2          -- Integer values of the number of faces for objects 1 and 2
                      respectively.

   *normals1
   *normals2       -- Arrays of doubles containing the co-ordinates of the
                      outward normal vectors to the faces of objects 1 and
                      2 respectively.  Normals are defined in the X, Y, Z
                      order with respect to a local coordinate frame of
                      reference. The object in this local coordinate frame
                      must be positioned such that it contains the origin.
                      The origin is also known as the seedpoint.

                      These normal arrays should be filled according to the
                      following syntax:

                         normals? = [x0, y0, z0, ... , xn, yn, zn]

                      Size of these arrays should be at least (3 * faces1)
                      and (3 * faces2) respectively.
   *xBars1
   *xBars2         -- Arrays of doubles containing the coordinates of points
                      on the faces of object 1 and 2 respectively.  Any point
                      on a face may be used.  Note that correspondence 
                      with the normal vector declared in "normals?" 
                      should be maintained.  Like the surface
                      normals, the defined coordinates are with respect to
                      the local coordinate frame with the origin being the
                      seedpoint.

                      These xBar arrays should be filled according to the
                      following syntax:

                         xBars? = [x0, y0, z0, ... , xn, yn, zn]

                      Size of these arrays should be at least (3 * faces1)
                      and (3 * faces2) respectively.

   *Gamma1
   *Gamma2         -- Arrays of doubles containing the co-ordinates of the
                      adjusted outward normal vectors to the faces of objects 1 and
                      2 respectively. These values are given as the output of the
                      routine 'setup_gamma'. Note that 'setup_gamma' is invoked only
                      once for each object. It uses normals and xBars defined above and
                      produces the adjusted outward normal as its output. Refer to 
                      documentations for 'setup_gamma' in 'interface.c' for details.
  
   *adjFaceCount1,
   *adjFaceCount2  -- Arrays of integers containing a count of the number
                      of adjacent faces to each face of objects 1 and 2
                      respectively.  Size of the array is equal to [faces1]
                      and [faces2] respectively.  For example,
                      adjFaceCount1[0] = 4 means that the first face of
                      object 1 has 4 adjacent faces.

   *adjSurf1,
   *adjSurf2       -- Arrays of integers containing the adjacent facial
                      information.  The size of this array should be
                      [faces1][adjCol1] and [faces2][adjCol2] respectively,
                      "adjCol?" is set to (1 + MAX(adjFaceCount?)).

                      NOTE: the face numbering MUST start from 0, not 1
                      and MUST be ordered either in a clockwise or
                      anti-clockwise manner.  In addition, each row
                      must be terminated by a -1.

                         Face ID : Adjacent faces list terminated by -1
                         -------   ------------------------------------
                            0    :      [1  2  5  8  -1 ]
                            1    :      [2  0  4  -1    ]
                            2    :      [0  1  3  -1    ]
                                               :
                                               :
                          faces1 :      [ ...           ]
                      
                      For example, face 0 has 4 adjacent faces whose ids are
                      1,2,5 and 8.
 
   
   adjCol1,
   adjCol2         -- Integer values that defines the size of the second
                      array subscript of "adjSurf1" and "adjSurf2".
   r1[3][3],
   r2[3][3]        -- 3 x 3 matrix of type double containing the rotational
                      transformation that will be applied to objects 1 and 2
                      respectively.  Both "r1" and "r2" specifies the transpose
                      of the rotation matrix to be used to rotate the objects.
                      The rotation matrix should be specified according to
                      the following notation:

                                      | R00  R10  R20 |
                                 r? = | R01  R11  R21 |
                                      | R02  R12  R22 |

                      where r?[0][0] = R00
                            r?[0][1] = R10
                            r?[0][3] = R20
                                  :
                            r?[2][2] = R22
   p1[3],
   p2[3]           -- 3 x 1 vector of type double defining the translation
                      to be applied to objects 1 and 2 respectively.

   jIndex[4][2]    -- Contains indices into "Gamma" array
                      of current active constraints.  This index is used
                      internally.
   increment       -- Indicates if incremental code should be used.  A
                      value of 1 enables incremental code.  A 0 disables
                      it.
   dInvReq         -- Indicates if the inverse of the C matrix is required.
                      The inverse is useful for getting gradient information
                      of the Growth Distance with respect to motion parameters.
                      If gradients are not needed, set this parameter to 0,
                      otherwise this parameter should be set to 1.  See the
                      appendix below for additional information on the gradient
                      returned in this parameter.
   *dWorkArray     -- A double workarray of size
                      (18 * (faces1 + faces2) * sizeof(double)).  This array
                      need not be initialized.  It is possible to share this
                      work array with multiple pairs of objects.  The
                      size of this work array should then be set to
                      (18 * MAX(faces1 + faces2) * sizeof(double)) where
                      MAX(faces1 + faces2) is the maximum of the total number
                      of faces of all object pairs.
   *iWorkArray     -- An integer workarray of size
                      (((6 * (faces1 + faces2)) + 9) * sizeof(int)).  NOTE:
                      it is important that this array be initialized to 0
                      the first time it is used.  Subsequently, this array
                      should not be touched.  Unlike dWorkArray, this work
                      array should not be shared among different paris of
                      object.

*/
int main(void)
{
    /* For this demonstration, we will define two objects, a cube
       and a diamond. */
    

    double *normals1, *xBars1;
    double *normals2, *xBars2;

    /* Adjacent face arrays. */
    int *adjSurf1, *adjFaceCount1, adjCol1;
    int *adjSurf2, *adjFaceCount2, adjCol2;

	/* Gamma Information.*/
	double *Gamma1, *Gamma2;

    /* Work arrays. */
    int *iWorkArray;
    double *dWorkArray;

    /* Storage space for indices needed by Growth Distance solver. */
    int jIndex[4][2];

    /* Number of faces information. */
    int planes1, planes2;

    /* Rotation matrix. */
    double r1[3][3], r2[3][3];

    /* Translation vector. */
    double p1[3], p2[3];

    /* Variable to hold the results of the computation. */
    double results[4];

    /* Specify if incremental code is to be used. */
    int increment;

    /*  -1
       D   matrix.  If "dInvReq" is set to 1, the Growth Distance
                               -1
       solver will return the D   matrix to the calling routine.
    */
    int dInvReq;
    double dInv[4][4];


    /* ====================== */
    /* Variables required for the computation of rotation matrices. */
    double c1, c2, c3, s1, s2, s3;
    double euler1[3], euler2[3];

    /* Local variables. */
    register int i;
    int *aPtr, count;
    double dT, *ptr;


    /* ========================================== */
    /* STEP 1: Definition of objects.             */
    /* ========================================== */
    /* In this demonstration routine, we will be using a cube and a
       diamond. */

    /* A cube has 6 faces, so initialize "planes1" to 6. */
    planes1 = 6;

    /* Each face of the cube has 4 adjacent faces.  Since each row
       of the "adjSurf" array needs to be terminated by a -1, we
       initialize "adjCol1" to 5. */
    adjCol1 = 5;

    /* Allocate memory for storage of "normals", "xBars" and the adjacent
       index "adjSurf". */
    i = planes1 * 3 * sizeof(double);
    normals1 = (double *)malloc(i);
    xBars1 = (double *)malloc(i);
    adjSurf1 = (int *)malloc(planes1 * adjCol1 * sizeof(int));
    adjFaceCount1 = (int *)malloc(planes1 * sizeof(int));

    /* Definition of surface normals for a unit cube. */
    ptr = normals1;
    *ptr++ = 1.0;
    *ptr++ = 0.0;
    *ptr++ = 0.0;
    *ptr++ = 0.0;
    *ptr++ = 1.0;
    *ptr++ = 0.0;
    *ptr++ = -1.0;
    *ptr++ = 0.0;
    *ptr++ = 0.0;
    *ptr++ = 0.0;
    *ptr++ = -1.0;
    *ptr++ = 0.0;
    *ptr++ = 0.0;
    *ptr++ = 0.0;
    *ptr++ = 1.0;
    *ptr++ = 0.0;
    *ptr++ = 0.0;
    *ptr = -1.0;

    /* Fill in corresponding "xBars" for the unit cube. 
       We are using the vertices of the cube. Any other
       points on the surface may be used.  */
    ptr = xBars1;
    *ptr++ = 0.5;
    *ptr++ = 0.5;
    *ptr++ = 0.5;
    *ptr++ = 0.5;
    *ptr++ = 0.5;
    *ptr++ = 0.5;
    *ptr++ = -0.5;
    *ptr++ = -0.5;
    *ptr++ = -0.5;
    *ptr++ = -0.5;
    *ptr++ = -0.5;
    *ptr++ = -0.5;
    *ptr++ = 0.5;
    *ptr++ = 0.5;
    *ptr++ = 0.5;
    *ptr++ = -0.5;
    *ptr++ = -0.5;
    *ptr = -0.5;

    /* Adjacent face list of cube. */
    aPtr = adjSurf1;
    *aPtr++ = 1;
    *aPtr++ = 4;
    *aPtr++ = 3;
    *aPtr++ = 5;
    *aPtr++ = -1;
    *aPtr++ = 0;
    *aPtr++ = 5;
    *aPtr++ = 2;
    *aPtr++ = 4;
    *aPtr++ = -1;
    *aPtr++ = 1;
    *aPtr++ = 5;
    *aPtr++ = 3;
    *aPtr++ = 4;
    *aPtr++ = -1;
    *aPtr++ = 0;
    *aPtr++ = 4;
    *aPtr++ = 2;
    *aPtr++ = 5;
    *aPtr++ = -1;
    *aPtr++ = 0;
    *aPtr++ = 1;
    *aPtr++ = 2;
    *aPtr++ = 3;
    *aPtr++ = -1;
    *aPtr++ = 0;
    *aPtr++ = 3;
    *aPtr++ = 2;
    *aPtr++ = 1;
    *aPtr = -1;

	/* For a cube, the number of adjacent faces to each face is fixed at
       4.  To generate the "adjFaceCount1" array, we can simply fill each
       element of the array with 4.  However, the following code presents
       a more generic approach that will work regardless of the object used.
       The key idea revolves around the fact that each row of the "adjSurf1"
       array is terminated by a -1. */
    for (i = 0; i < planes1; i++)
    {
        /* Point to start of adjacent face list. */
        aPtr = adjSurf1 + (i * adjCol1);

        count = 0;
        while (*aPtr++ != -1)
            ++count;

        /* Update information. */
        adjFaceCount1[i] = count;
    }

    /* ====================== */
    /* A diamond has 8 faces, so initialize "planes2" to 8. */
    planes2 = 8;

    /* Each face of the diamond has 3 adjacent faces.  Since each row
       of the "adjSurf" array needs to be terminated by a -1, we
       initialize "adjCol2" to 4. */
    adjCol2 = 4;

    /* Allocate memory for storage of "normals", "xBars" and
       "adjSurf". */
    i = planes2 * 3 * sizeof(double);
    normals2 = (double *)malloc(i);
    xBars2 = (double *)malloc(i);
    adjSurf2 = (int *)malloc(planes2 * adjCol2 * sizeof(int));
    adjFaceCount2 = (int *)malloc(planes2 * sizeof(int));

    /* Definition of surface normals for a unit diamond. */
    ptr = normals2;
    *ptr++ = 1.0;
    *ptr++ = 1.0;
    *ptr++ = -1.0;
    *ptr++ = -1.0;
    *ptr++ = 1.0;
    *ptr++ = -1.0;
    *ptr++ = -1.0;
    *ptr++ = 1.0;
    *ptr++ = 1.0;
    *ptr++ = 1.0;
    *ptr++ = 1.0;
    *ptr++ = 1.0;
    *ptr++ = 1.0;
    *ptr++ = -1.0;
    *ptr++ = 1.0;
    *ptr++ = 1.0;
    *ptr++ = -1.0;
    *ptr++ = -1.0;
    *ptr++ = -1.0;
    *ptr++ = -1.0;
    *ptr++ = -1.0;
    *ptr++ = -1.0;
    *ptr++ = -1.0;
    *ptr = 1.0;

    /* Fill in corresponding "xBars" for the unit diamond. */
    ptr = xBars2;
    *ptr++ = 0.0;
    *ptr++ = 1.0;
    *ptr++ = 0.0;
    *ptr++ = 0.0;
    *ptr++ = 1.0;
    *ptr++ = 0.0;
    *ptr++ = 0.0;
    *ptr++ = 1.0;
    *ptr++ = 0.0;
    *ptr++ = 0.0;
    *ptr++ = 1.0;
    *ptr++ = 0;
    *ptr++ = 0.0;
    *ptr++ = -1.0;
    *ptr++ = 0.0;
    *ptr++ = 0.0;
    *ptr++ = -1.0;
    *ptr++ = 0.0;
    *ptr++ = 0.0;
    *ptr++ = -1.0;
    *ptr++ = 0.0;
    *ptr++ = 0.0;
    *ptr++ = -1.0;
    *ptr = 0;

    /* Adjacent face list of diamond. */
    aPtr = adjSurf2;
    *aPtr++ = 1;
    *aPtr++ = 3;
    *aPtr++ = 5;
    *aPtr++ = -1;
    *aPtr++ = 0;
    *aPtr++ = 6;
    *aPtr++ = 2;
    *aPtr++ = -1;
    *aPtr++ = 1;
    *aPtr++ = 7;
    *aPtr++ = 3;
    *aPtr++ = -1;
    *aPtr++ = 0;
    *aPtr++ = 2;
    *aPtr++ = 4;
    *aPtr++ = -1;
    *aPtr++ = 3;
    *aPtr++ = 7;
    *aPtr++ = 5;
    *aPtr++ = -1;
    *aPtr++ = 0;
    *aPtr++ = 4;
    *aPtr++ = 6;
    *aPtr++ = -1;
    *aPtr++ = 1;
    *aPtr++ = 5;
    *aPtr++ = 7;
    *aPtr++ = -1;
    *aPtr++ = 2;
    *aPtr++ = 6;
    *aPtr++ = 4;
    *aPtr = -1;

    /* An approach similar to that of the cube may be used to fill up
       the "adjFactCount2" array for the diamond. */
    for (i = 0; i < planes2; i++)
    {
        /* Point to start of adjacent face list. */
        aPtr = adjSurf2 + (i * adjCol2);

        count = 0;
        while (*aPtr++ != -1)
            ++count;

        /* Update information. */
        adjFaceCount2[i] = count;
    }

    /* ====================== */
    /* Finally, we have to allocate sufficient space for the work arrays
       needed by the Growth Distance solver. */
    dWorkArray = (double *)malloc(18 * (planes1 + planes2) * sizeof(double));
    iWorkArray = (int *)malloc(((6 * (planes1 + planes2)) + 9) * sizeof(int));

    /* We have to initialize the "iWorkArray" to 0.  NOTE: we need to do this
       only once.  Subsequently, when the Growth Distance solver is running,
       we should not modify this array. */
    for (i = 0; i < ((6 * (planes1 + planes2)) + 1); i++)
        iWorkArray[i] = 0;


    /* ========================================== */
    /* STEP 2:  Preparing the data.               */
    /* ========================================== */
    
    /* Setup gamma information for each object.  This step is neccesary
	   for solver function to run correctly.*/
    setup_gamma(normals1, xBars1, planes1, &Gamma1 );
    setup_gamma(normals2, xBars2, planes2, &Gamma2 );
	
    /* Form the rotation matrix for object 1.  Here we use the Euler
       angle rotational matrix defined by

                  |      c2.c3           -c2.s3        s2   |
              r = | s1.s2.c3+c1.s3  -s1.s2.s3+c1.c3  -s1.c2 |
                  |-c1.s2.c3+s1.s3   c1.s2.s3+s1.c3   c1.c2 |

       We begin by defining the rotation and translation parameters,
       "euler1" and "p1" for object 1 first.
    */
    euler1[0] = 10.0 * PIOVER180;
    euler1[1] = 45.0 * PIOVER180;
    euler1[2] = -33.0 * PIOVER180;
    p1[0] = 0.0;
    p1[1] = 0.0;
    p1[2] = 0.0;

    /* Compute rotation matrix for object 1. */
    c1 = cos(euler1[0]);
    c2 = cos(euler1[1]);
    c3 = cos(euler1[2]);
    s1 = sin(euler1[0]);
    s2 = sin(euler1[1]);
    s3 = sin(euler1[2]);
    r1[0][0] = c2 * c3;
    r1[0][1] = -(c2 * s3);
    r1[0][2] = s2;
    r1[1][0] = s1 * s2 * c3 + c1 * s3;
    r1[1][1] = c1 * c3 - s1 * s2 * s3;
    r1[1][2] = -(s1 * c2);
    r1[2][0] = s1 * s3 - c1 * s2 * c3;
    r1[2][1] = c1 * s2 * s3 + s1 * c3;
    r1[2][2] = c1 * c2;

    /* Define the rotation and translation parameters, "euler2" and "p2"
       for object 2. */
    euler2[0] = 2.0 * PIOVER180;
    euler2[1] = 5.0 * PIOVER180;
    euler2[2] = -80.0 * PIOVER180;
    p2[0] = 10.0;
    p2[1] = 0.0;
    p2[2] = 0.0;

    /* Compute rotation matrix for object 2. */
    c1 = cos(euler2[0]);
    c2 = cos(euler2[1]);
    c3 = cos(euler2[2]);
    s1 = sin(euler2[0]);
    s2 = sin(euler2[1]);
    s3 = sin(euler2[2]);
    r2[0][0] = c2 * c3;
    r2[0][1] = -(c2 * s3);
    r2[0][2] = s2;
    r2[1][0] = s1 * s2 * c3 + c1 * s3;
    r2[1][1] = c1 * c3 - s1 * s2 * s3;
    r2[1][2] = -(s1 * c2);
    r2[2][0] = s1 * s3 - c1 * s2 * c3;
    r2[2][1] = c1 * s2 * s3 + s1 * c3;
    r2[2][2] = c1 * c2;

    /*                     -1
       If you require the D   matrix, set "dInvReq" to 1.  For this
       demonstration, we set dInvReg to 1, or set dInvReg to zero.*/
    dInvReq = 1;

    /* The first time we call the Growth Distance solver, we should not
       use incremental code, so we set "increment" to 0. */
    increment = 0;

    /* ====================== */
    /* Print header message. */
    printf("\n    X\t\t    Y\t\t    Z\t\t Sigma");
    printf("\n   ===\t\t   ===\t\t   ===\t\t=======\n");


    /* ============================================ */
    /* STEP 3:  Calling the Growth Distance Solver. */
    /* ============================================ */
    dT = 0.5;
    for (i = 0; i < 10; i++)
    {
        /* Startup Growth Distance Solver. */
        gd( Gamma1, Gamma2,adjSurf1, adjSurf2, adjFaceCount1, adjFaceCount2,
			adjCol1, adjCol2, r1, r2, p1, p2, jIndex, results, dInv, dInvReq, increment,
			planes1, planes2, dWorkArray, iWorkArray);

        /* Print out resutls. */
        printf("%f\t%f\t%f\t%f\n", results[0], results[1], results[2],
								   results[3]);

        /* Move object 2 a little. */
        p2[0] += dT;
        p2[1] += dT;
        p2[2] += dT;

        /* Let us set the "increment" flag so that the next time we call the
           Growth Distance solver, we make use of the faster incremental
           code. */
        increment = 1;
    }

    return 0;
}
