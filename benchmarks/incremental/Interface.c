#include "Gamma.h"
#include <malloc.h>
#include "Solver4c.h"


/* ================================================================ */
/* This is the entry point into this Growth Distance main function 'gd'.  It
   initializes all variables and prepares the necessary data required
   by the main function.

   INPUT PARAMETERS TO 'GD' FUNCTION
   =================================

   faces1,
   faces2          -- Integer values of the number of faces for objects 1 and 2
                      respectively.
   
   *Gamma1
   *Gamma2         -- Arrays of doubles containing the co-ordinates of the
                      adjusted outward normal vectors to the faces of objects 1 and
                      2 respectively. These values should be the output of the
                      routine setup_gamma. Note that setup_gamma is invoked only
                      once for each object. The main routine 'gd' for growth distance
                      computation calls to 'gd' use Gamma1 and Gamma2 directly. See
                      the following documentations for setup_gamma for details.
  
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

                      NOTE: the face numbering must start from 0, not 1
                      and must be ordered either in a clockwise or
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


   OUTPUT PARAMETERS FROM 'GD' FUNCTION
   ====================================
   result[4]       -- The first 3 x 1 vector (result[0] .. result[2])
                      contains the coordinates of the point of contact
                      of the objects corresponding to the sigma value
                      stored in "result[3]".
   dInv[4][4]      -- If "dInvReq" is set to 1, this 4 x 4 double array
                      returns the inverse of the C matrix.  See the appendix
                      below on how dInv can be used for gradient computations.


   RETURN VALUE OF 'GD' FUNCTION
   ==============================
   The return values of this function may be interpreted as follows:
        1:    Optimality has been reached.  The solution is returned in
              the variable "result".
        2:    The two seed points are coincident.  Solution has been
              reached.
        3:    Divide by zero error.  Two possible causes are:
               i. Possible error in definition of input data.
              ii. Face normals of adjacent faces are almost parallel.
                  This could mean that representation of the object
                  is too fine, for example, representing a unit sphere
                  with 10,000 faces.

========
APPENDIX
========
The solution for Growth Distance is formulated as an LP:

                      Minimize c'z

                s.t
	                    Ez <= b

where c = [0 0 0 1]', z = [x' sigma]', and E contains the coefficients of
x and sigma corresponding to the facial constraints of the two objects.  On
output, the optimum x and sigma are defined by 4 active constraints.  The
collection of the coefficients of these 4 constraints is C such that
                          _   _
                         Cz = b
      _    _   _____
where z = [x'  sigma]' is the optimum solution of the LP.  The inverse of C
provides gradient information with respect to configuration parameters theta
(hereafter theta will be represented by the symbol t) defining objects A and B.
The expression for the gradient with respect to theta is 
                                _
              dz        -1     db       dC    _
              --(t ) = C  (t )[--(t ) - --(t )z(t )]
              dt  0         0  dt  0    dt  0    0

                                               -1
Choosing dInvReq as 1 will provide the matrix C  (t ).  The expressions
 _                                                 0
db/dt and dC/dt follows from the kinematic relationship of b and C on t.  For
details, please refer to "Growth Distances: New Measures for Object Separation
and Penetration" by Chong-Jin Ong and Elmer G. Gilbert, IEEE Transactions on
Robotics and Automation, Vol 12, No. 6, December 1996, pg 898.
*/
int gd(double * Gamma1, double * Gamma2,int *adjSurf1, int *adjSurf2,
       int *adjFaceCount1, int *adjFaceCount2, int adjCol1,
       int adjCol2, double r1[][3], double r2[][3], double p1[],
       double p2[], int jIndex[][2], double result[],
	   double dInv[][4], int dInvReq, int increment, int faces1,
	   int faces2, double *dWorkArray, int *iWorkArray)
{
	int *nA, *nB;
	int *iIndex1, *iIndex2;

	/* Initialize pointers into work array. */
	nA = iWorkArray;
	nB = nA + 1;
	iIndex1 = nB + 1;
	iIndex2 = iIndex1 + 3;

    /* Initialize pointers to gamma arrays. */
    initGammaArrays(faces1, faces2, dWorkArray, iIndex2 + 3);

    /* This is the main call to the solver routine. */
    return solver(Gamma1, Gamma2, adjSurf1,
                  adjSurf2, adjFaceCount1, adjFaceCount2, adjCol1,
                  adjCol2, r1, r2, p1, p2, iIndex1, iIndex2, nA,
                  nB, jIndex, result, dInv, dInvReq, increment);
}

/* ================================================================ */
/* This function setup_gamma generates gamma information as required by 'gd'. 

   INPUT PARAMETERS
   ================
   faces           -- Integer values of the number of faces for this object.
   
   *normals        -- Array of doubles containing the co-ordinates of the
                      outward normal vectors to the faces of this object.
					  Normals are defined in the X, Y, Z
                      order with respect to a local coordinate frame of
                      reference, the origin, which is know as the seed point.

                      This normal array should be filled according to the
					  following syntax:

                         normals? = [x0, y0, z0, ... , xn, yn, zn]

                      Size of these arrays should be at least (3 * faces)
                      and (3 * faces) respectively.

                        
   *xBars          -- Array of doubles containing the coordinates of points
                      on the faces of this object.  Note that
                      correspondence with the normal vector declared in
                      normals should be maintained.  

                      This xBar array should be filled according to the
					  following syntax:

                           xBars? = [x0, y0, z0, ... , xn, yn, zn]

                      Size of these arrays should be at least (3 * faces)
                      and (3 * faces) respectively.

   OUTPUT PARAMETERS
   =================
   
   *Gamma          -- Array of doubles containing the co-ordinates of the
                      adjusted outward normal vectors to the faces of this object.
                      This size of this array should be the same size as the
                      normals parameter
                             
*/

void setup_gamma(double *normals, double *xBars, int faces, double **Gamma)
{
	
	int i ;
	double bottom ;
	double *nor_Ptr, * xBar_Ptr, *Gamma_Ptr;
      double nor_x, nor_y, nor_z;
	double xBar_x, xBar_y, xBar_z;

	
	/* Compute memory size needed for Gamma.*/
    i = faces * 3 * sizeof(double);
    *Gamma = (double *)malloc(i);
	
	/*give the three pointer original address. */
	nor_Ptr   = normals ; 
	xBar_Ptr  = xBars ; 
	Gamma_Ptr = *Gamma ;


    for(i = 0; i< faces; i++)
	{   
        nor_x = *nor_Ptr++;
        nor_y = *nor_Ptr++;
        nor_z = *nor_Ptr++;

        xBar_x = *xBar_Ptr++;
		xBar_y = *xBar_Ptr++;
		xBar_z = *xBar_Ptr++;
        

        /*                       1
           "Bottom" is equal to ---- .
                                  -
                                n.x
        */
        bottom = 1.0 / ( nor_x * xBar_x + nor_y * xBar_y + nor_z * xBar_z);

        /*            
                     n    
           Compute  ---- .
                      -
                    n.x
        */
        *Gamma_Ptr++ =  nor_x * bottom;
		*Gamma_Ptr++ =  nor_y * bottom;    
		*Gamma_Ptr++ =  nor_z * bottom;
    
        
    }

}


