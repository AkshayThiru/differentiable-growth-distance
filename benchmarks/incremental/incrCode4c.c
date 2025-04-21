#include <stdio.h>
#include <math.h>
#include "Gamma.h"
#include "Tolerance.h"
#include "incrCode4c.h"
#include "Solver4c.h"


/* External variables. */
extern double gammaVals[4][3];


/* Local function prototypes. */
__inline void getFeasiblePoint(void (*getGamma)(),double *Gamma, double r[][3],
                               double xMinusP[], int iIndex[],
                               int violatedFace, double result[]);
__inline void recover(void (*getGammaFailed)(), double *GammaFailed,
					  int *adjSurfFailed, int adjColFailed, double rFailed[][3],
                      double xMinusPFailed[], int iIndexFailed[],
                      int *nFailed, int violatedFace, int failedObject,
                      void (*getGamma)(), double *Gamma, int *adjSurf, int adjCol,
                      double r[][3], double xMinusP[], int iIndex[],
                      int *n, int jIndex[][2], double result[]);
__inline int linearSolver3x3(double A[][3], double b[], double x[]);



/* ================================================================ */
/* This function implements incremental coding for improved
   efficiency in obtaining an optimal solution.


   Input Parameters
   ================
   *Gamma2
   *Gamma1  -- Arrays of doubles containing the co-ordinates of the
                      outward original gamma vectors to the faces of objects 1 and
                      2 respectively.  origin_Gamma are defined in the X, Y, Z
                      order with respect to a local coordinate frame of
                      reference, the origin, which is know as the seed point.

                         origin_Gamma? = [x0, y0, z0, ... , xn, yn, zn]

                      Size of these arrays should be at least (3 * faces1)
                      and (3 * faces2) respectively.

					  When you write your own program, you need call setup_origin_gamma( )
					  to create origin_Gamma array for each object with respect to
					  the values of normals and xBars you input.

   *adjSurf1,
   *adjSurf2     -- Arrays containing the adjacent surface information.
                    The size of this array should be [planes1][adjCol1]
                    and [planes2][adjCol2] respectively.  Note that
                    surface numbering must start from 0, not 1.
   adjCol1,
   adjCol2       -- Defines the size of the second array subscript of
                    "adjSurf1" and "adjSurf2".
   r1[3][3],
   r2[3][3]      -- Specifies the transpose of the rotation matrix to
                    be used to rotate the objects.  The rotation matrix
                    should be specified according to the notation

                                    | R11  R21  R31 |
                               r? = | R12  R22  R32 |
                                    | R13  R23  R33 |

                    where r?[0][0] = R11
                          r?[0][1] = R21
                          r?[0][3] = R31
                                :
                          r?[2][2] = R33

                    Use these matrices to rotate the objects with
                    respect to their seed points.
   p1[3],
   p2[3]         -- These 3 x 1 vectors specifices the translation
                    vectors associated with each object respectively.
                    Use these vectors to position the seed points
                    of the objects with respect to the global frame
                    of reference.
   iIndex1[3],
   iIndex2[3]    -- Index of active faces for objects 1 and 2.
   *nA, *nB      -- Count of number of elements in "iIndex1" and
                    "iIndex2" respectively.
   jIndex[4][2]  -- Contains indices into "normals" and "xBars" arrays
                    of current active constraints.


   Output Parameter
   ================
   result[4]     -- The first 3 x 1 vector (result[0] .. result[2])
                    contains the x coordinate of the point of contact
                    of the objects corresponding to the sigma value
                    stored in "result[3]".


   Return Value
   ============
   The return values of this function may be interpreted as follows:
        0:    Feasibility is maintained.
        1:    Feasibility is not maintained and Delta Angle is smaller than sth,
		                 Recovery attempted.
        2:    Seed points are coincident.  Solution has been reached.
        3:    Feasibility is not maintained and Delta Angle is greater than sth,
		              or error occurring, Do restarting.
*/
int incrSolver(double *Gamma1, double *Gamma2,
			   int *adjSurf1, int *adjSurf2,
               int *adjFaceCount1, int *adjFaceCount2,
               int adjCol1, int adjCol2, double r1[][3],
               double r2[][3], double p1[], double p2[],
               int iIndex1[], int iIndex2[], int *nA, int *nB,
               int jIndex[][2], double result[])
{
   
	/* Pointer to adjSurf arrays.*/
    int *adjPtr;

    /* Variables used to find the new x' and growth factor sigma'. */
    double A[3][3], b[3], gammaP[4];

    /* Variables to keep track of face violations. */
    int violatedFace1, violatedFace2, violatedObject, noOfViolatedFace;
    double violatedFaceValue;

    /* Miscellaneous variables. */
    int i,  active;
    double value, reference;
    double xMinusPA[3], origin_xMinusPA[3], xMinusPB[3], origin_xMinusPB[3];
	double returnVal[3], *Gamma_Ptr;

	 /* ====================== */
#ifdef _CHECK_SEED
    double dP[3];

    /* Get direction vector joining the two seed points. */
    dP[0] = p2[0] - p1[0];
    dP[1] = p2[1] - p1[1];
    dP[2] = p2[2] - p1[2];

	/* Check if the seed points are coincident. */
    if (((dP[0] > -TOLERANCE) && (dP[0] < TOLERANCE)) &&
        ((dP[1] > -TOLERANCE) && (dP[1] < TOLERANCE)) &&
        ((dP[2] > -TOLERANCE) && (dP[2] < TOLERANCE)))
    {
        /* The two seed points conincide.  In theory, the growth
           factor "sigma" is infinite.  To cater for this, we return
           a zero for the growth function "sigma" as the solution. */
        result[0] = p1[0];
        result[1] = p1[1];
        result[2] = p1[2];
        result[3] = 0.0;

        return 2;
    }

#endif

    /* Get gamma values for the active constraints. */
    active = *nA + *nB;

    switch (active)
    {
        case 0:
            /* If the number of active constraints is zero, and assuming
               the user did initialize the variables "nA" and "nB" once
               at the start of the driver routine, the user has called
               this incremental code without a valid set of constraints.
               So we will restart the algorithm from the top. */
        case 2:
        case 3:
            /* With only 2 or 3 active constraints, we cannot compute a
               value for x' and gamma'.  So just restart from scratch
               but with a better guess of which face to start the
               search. */
			{	return 3;}

        case 4:
            /* We have 4 active constraints. */
            for (i = 0; i < 4; i++)
            {
                if (jIndex[i][1] == 1)
                    getGamma1(Gamma1, r1, jIndex[i][0],
                              gammaVals[i]);
                else
                    getGamma2(Gamma2, r2, jIndex[i][0],
                              gammaVals[i]);
            }

            /* Need to find the new x' and growth factor sigma'.
                                    T         T
               First compute (gammaA  - gammaB ).
            */
            A[0][0] = gammaVals[0][0] - gammaVals[1][0];
            A[0][1] = gammaVals[0][1] - gammaVals[1][1];
            A[0][2] = gammaVals[0][2] - gammaVals[1][2];

            /*                T         T
               Compute (gammaA  - gammaC ).
            */
            A[1][0] = gammaVals[0][0] - gammaVals[2][0];
            A[1][1] = gammaVals[0][1] - gammaVals[2][1];
            A[1][2] = gammaVals[0][2] - gammaVals[2][2];

            /*                T         T
               Compute (gammaA  - gammaD ).
            */
            A[2][0] = gammaVals[0][0] - gammaVals[3][0];
            A[2][1] = gammaVals[0][1] - gammaVals[3][1];
            A[2][2] = gammaVals[0][2] - gammaVals[3][2];

            /*                T
               Compute (gamma? .P ) for all 4 constraints
                                 ?
            */
            for (i = 0; i < 4; i++)
            {
                if (jIndex[i][1] == 1)
                {
                    /* Gamma value from object 1. */
                    gammaP[i] = getGamma1P1(Gamma1, r1, p1,
                                            jIndex[i][0]);
                }
                else
                {
                    /* Gamma value from object 2. */
                    gammaP[i] = getGamma2P2(Gamma2, r2, p2,
                                            jIndex[i][0]);
                }
            }

            /* Ready "b" matrix.  First compute
                      T            T
               (gammaA .P  - gammaB .P  )
                         A            B
            */
            b[0] = gammaP[0] - gammaP[1];


            /*                T            T
               Compute (gammaA .P  - gammaC .P  )
                                 A            C
            */
            b[1] = gammaP[0] - gammaP[2];


            /*                T            T
               Compute (gammaA .P  - gammaD .P  )
                                 A            D
            */
            b[2] = gammaP[0] - gammaP[3];

            /* Solve for x'. */
            if (linearSolver3x3(A, b, result) == 1)
                /* Have to restart. */
			{ return 3;}

            /*                          T
               Solve for sigma' = gammaA (x' - P ).
                                                A
            */
            result[3] = gammaVals[0][0] * result[0] +
                        gammaVals[0][1] * result[1] +
                        gammaVals[0][2] * result[2] - gammaP[0];

            /* "sigma" value cannot be negative. */
            if (result[3] < -TOLERANCE)
                /* We are way out already, have to restart. */
            {   return 3;}

    }

    /* Store (x' - P ).
                    A
    */
    xMinusPA[0] = result[0] - p1[0];
    xMinusPA[1] = result[1] - p1[1];
    xMinusPA[2] = result[2] - p1[2];

    /* Store (x' - P ).
                    B
    */
    xMinusPB[0] = result[0] - p2[0];
    xMinusPB[1] = result[1] - p2[1];
    xMinusPB[2] = result[2] - p2[2];

    /* ========================================== */
    /* Set reference point for comparison. */
    reference = result[3] ;

    /* Initialize variables. */
    noOfViolatedFace = 0;

    /* Check for feasibility. */
    if (*nA == 3)
    {
        /* We have a vertex-face contact.  We only need to check
           face of object 2 for feasibility. */
        adjPtr = adjSurf2 + iIndex2[0] * adjCol2;

	    /*rotate the xMinusP to the original coordinate .*/
        get_origin_vector(xMinusPB, origin_xMinusPB, r2); 
		do
		{
            /*                                       -1
            Compute "bottom" = Gamma .dx  . */
            Gamma_Ptr = Gamma2 + 3*(*adjPtr);
		    returnVal[0] = * Gamma_Ptr ++;
		    returnVal[1] = * Gamma_Ptr ++;
		    returnVal[2] = * Gamma_Ptr ;
                
			value = returnVal[0] * origin_xMinusPB[0] +
                        returnVal[1] * origin_xMinusPB[1] +
                        returnVal[2] * origin_xMinusPB[2];

            if (value > reference)
			{   
				 
				    /* New x' and sigma' is not feasible, so attempt
                    a recovery.  */
                    recover(getGamma2, Gamma2,adjSurf2,
                        adjCol2, r2, xMinusPB, iIndex2,
                        nB, *adjPtr, 2, getGamma1, 
                        Gamma1, adjSurf1, adjCol1, r1, xMinusPA,
                        iIndex1, nA, jIndex, result);

                    /* Return appropriate code. */
                    return 1;
			}

            /* Next adjacent surface. */
            adjPtr++;

		} while (*adjPtr != -1);
    }
    else if (*nB == 3)
    {
        /* We have a face-vertex contact.  We only need to check
           object 1 for feasibility. */
        adjPtr = adjSurf1 + iIndex1[0] * adjCol1;

       /*rotate the xMinusP to the original coordinate .*/
       get_origin_vector(xMinusPA, origin_xMinusPA, r1); 
	   do
	 {
	         
		     /*                                       -1
             Compute "bottom" = Gamma .dx  . */
             Gamma_Ptr = Gamma1 + 3*(*adjPtr);
		     returnVal[0] = * Gamma_Ptr ++;
		     returnVal[1] = * Gamma_Ptr ++;
		     returnVal[2] = * Gamma_Ptr ;

             value = returnVal[0] * origin_xMinusPA[0] +
                           returnVal[1] * origin_xMinusPA[1] +
                           returnVal[2] * origin_xMinusPA[2];

            if (value > reference)
			{    
					   
					/* New x' and sigma' is not feasible, so attempt
                          a recovery. */
                            recover(getGamma1, Gamma1,adjSurf1,
                                        adjCol1, r1, xMinusPA, iIndex1,
                                        nA, *adjPtr, 1, getGamma2, 
                                        Gamma2,adjSurf2, adjCol2, r2, xMinusPB,
                                        iIndex2 , nB, jIndex, result);

                               /* Return appropriate code. */
                               return 1;
				 
			}

          /* Next adjacent surface. */
          adjPtr++;

		} while (*adjPtr != -1);
		
    }
    else
    {
        /* Initialize variables. */
        violatedFaceValue = 1e10;

        if (adjFaceCount1[iIndex1[0]] <= adjFaceCount1[iIndex1[1]])
        {
            /* We have edge-edge contact.  We have no choice but to
               check all neighouring faces of one face of each object.
			   (We choose one face from two active faces which has 
			   fewer neighbouring faces.) 
			   We start by computing index into adjacent list of current
			   surface. This current surface is given by the "iIndex1".*/
            adjPtr = adjSurf1 + iIndex1[0] * adjCol1;


            /* This version does not take into account that some
               adjacent surfaces may overlap.  It just make sure
               that none of the faces tested belongs to the active
               set. */
            /*rotate the xMinusP to the original coordinate .*/
            get_origin_vector(xMinusPA, origin_xMinusPA, r1); 

			do
            {
                if (*adjPtr != iIndex1[1])
                {
                    /*                                       -1
                          Compute "bottom" = Gamma .dx  . */
                    Gamma_Ptr = Gamma1 + 3*(*adjPtr);
			        returnVal[0] = * Gamma_Ptr ++;
			        returnVal[1] = * Gamma_Ptr ++;
			        returnVal[2] = * Gamma_Ptr ;

                    value = returnVal[0] * origin_xMinusPA[0] +
                            returnVal[1] * origin_xMinusPA[1] +
                            returnVal[2] * origin_xMinusPA[2];

                    if (value > reference)
                    {
                        /* New x' and sigma' is not feasible. */
                        violatedFaceValue = value;
                        violatedFace1 = *adjPtr;
                        violatedObject = 1;
                        ++noOfViolatedFace;

                        break;
                    }
                }

                /* Next adjacent surface. */
                adjPtr++;

            } while (*adjPtr != -1);
        }
        else
        {
             /* We have edge-edge contact.  We have no choice but to
               check all neighouring faces of one face of each object.
			   (We choose one face from two active faces which has 
			   fewer neighbouring faces.) 
			   We start by computing index into adjacent list of current
			   surface. This current surface is given by the "iIndex1".*/
            adjPtr = adjSurf1 + iIndex1[1] * adjCol1;

            /* This version does not take into account that some
               adjacent surfaces may overlap.  It just make sure
               that none of the faces tested belongs to the active
               set. */

			/*rotate the xMinusP to the original coordinate .*/
            get_origin_vector(xMinusPA, origin_xMinusPA, r1); 
            do
            {
                if (*adjPtr != iIndex1[0])
                {

		            /*                                       -1
                          Compute "bottom" = Gamma .dx  . */
                    Gamma_Ptr = Gamma1 + 3*(*adjPtr);
			        returnVal[0] = * Gamma_Ptr ++;
			        returnVal[1] = * Gamma_Ptr ++;
			        returnVal[2] = * Gamma_Ptr ;
                    
					value = returnVal[0] * origin_xMinusPA[0] +
                            returnVal[1] * origin_xMinusPA[1] +
                            returnVal[2] * origin_xMinusPA[2];

                    if (value > reference)
                    {
                        /* New x' and sigma' is not feasible. */
                        violatedFaceValue = value;
                        violatedFace1 = *adjPtr;
                        violatedObject = 1;
                        ++noOfViolatedFace;

                        break;
                    }
                }
                
                /* Next adjacent surface. */
                adjPtr++;


            } while (*adjPtr != -1);
        }

        if (adjFaceCount2[iIndex2[0]] <= adjFaceCount2[iIndex2[1]])
        {
            /* ============================ */
            /* Check the other surface from object 2 for feasibility. */
			/* We have edge-edge contact.  We have no choice but to
               check all neighouring faces of one face of each object.
			   (We choose one face from two active faces which has 
			   fewer neighbouring faces.) 
			   We start by computing index into adjacent list of current
			   surface. This current surface is given by the "iIndex1".*/
            adjPtr = adjSurf2 + iIndex2[0] * adjCol2;


            /* This version does not take into account that some
               adjacent surfaces may overlap.  It just make sure
               that none of the faces tested belongs to the active
               set. */
	        /*rotate the xMinusP to the original coordinate .*/
            get_origin_vector(xMinusPB, origin_xMinusPB, r2); 

            do
            {
                if (*adjPtr != iIndex2[1])
                {
              
		           /*                                       -1
                         Compute "bottom" = Gamma .dx  . */
                   Gamma_Ptr = Gamma2 + 3*(*adjPtr);
			       returnVal[0] = * Gamma_Ptr ++;
			       returnVal[1] = * Gamma_Ptr ++;
			       returnVal[2] = * Gamma_Ptr ;
	
                    value = returnVal[0] * origin_xMinusPB[0] +
                            returnVal[1] * origin_xMinusPB[1] +
                            returnVal[2] * origin_xMinusPB[2];

                    if (value > reference)
                    {
                        /* New x' and sigma' is not feasible. */
                        ++noOfViolatedFace;
                        violatedFace2 = *adjPtr;
                        if (value < violatedFaceValue)
                            violatedObject = 2;

                        break;
                    }
                }

                /* Next adjacent surface. */
                adjPtr++;

            } while (*adjPtr != -1);
        }
        else
        {
            /* ============================ */
            /* Check the other surface from object 2 for feasibility. */
			/* We have edge-edge contact.  We have no choice but to
               check all neighouring faces of one face of each object.
			   (We choose one face from two active faces which has 
			   fewer neighbouring faces.) 
			   We start by computing index into adjacent list of current
			   surface. This current surface is given by the "iIndex1".*/
            adjPtr = adjSurf2 + iIndex2[1] * adjCol2;

            /* This version does not take into account that some
               adjacent surfaces may overlap.  It just make sure
               that none of the faces tested belongs to the active
               set. */

	        /*rotate the xMinusP to the original coordinate .*/
            get_origin_vector(xMinusPB, origin_xMinusPB, r2); 

            do
            {
                if (*adjPtr != iIndex2[0])
                {

		            /*                                       -1
                          Compute "bottom" = Gamma .dx  . */
                    Gamma_Ptr = Gamma2 + 3*(*adjPtr);
			        returnVal[0] = * Gamma_Ptr ++;
			        returnVal[1] = * Gamma_Ptr ++;
			        returnVal[2] = * Gamma_Ptr ;

                    value = returnVal[0] * origin_xMinusPB[0] +
                            returnVal[1] * origin_xMinusPB[1] +
                            returnVal[2] * origin_xMinusPB[2];

                    if (value > reference)
                    {
                        /* New x' and sigma' is not feasible. */
                        ++noOfViolatedFace;
                        violatedFace2 = *adjPtr;
                        if (value < violatedFaceValue)
                            violatedObject = 2;

                        break;
                    }
                }

                /* Next adjacent surface. */
                adjPtr++;

            } while (*adjPtr != -1);
        }
    }

    if (noOfViolatedFace == 0)
        /* The new computed values of x' and sigma' is feasible.
           So return with no error. */
        return 0;
    else 
	{if (noOfViolatedFace == 1)
	  {
        if (violatedObject == 1)
        {
            /* New x' and sigma' is not feasible, so attempt
               a recovery. */
            recover(getGamma1, Gamma1,adjSurf1,
                    adjCol1, r1, xMinusPA, iIndex1,
                    nA, violatedFace1, 1, getGamma2, 
                    Gamma2, adjSurf2, adjCol2, r2, xMinusPB,
                    iIndex2, nB, jIndex, result);
        }
        else
        {
            /* New x' and sigma' is not feasible, so attempt
               a recovery. */
            recover(getGamma2, Gamma2,adjSurf2,
                    adjCol2, r2, xMinusPB, iIndex2,
                    nB, violatedFace2, 2, getGamma1,
                    Gamma1, adjSurf1, adjCol1, r1, xMinusPA,
                    iIndex1, nA, jIndex, result);
        }
	  }
     else
	 {
        if (violatedObject == 1)
        {
            /* New x' and sigma' is not feasible, so attempt
               a recovery. */
            getFeasiblePoint(getGamma1, Gamma1,r1,
                             xMinusPA, iIndex1, violatedFace1,
                             result);

            /* Update (x - P ).
                            A
            */
            xMinusPA[0] = result[0] - p1[0];
            xMinusPA[1] = result[1] - p1[1];
            xMinusPA[2] = result[2] - p1[2];

            /* Update (x - P ).
                            B
            */
            xMinusPB[0] = result[0] - p2[0];
            xMinusPB[1] = result[1] - p2[1];
            xMinusPB[2] = result[2] - p2[2];

            /* Recover from infeasibility. */
            recover(getGamma2, Gamma2, adjSurf2,
                    adjCol2, r2, xMinusPB, iIndex2,
                    nB, violatedFace2, 2, getGamma1,
                    Gamma1,adjSurf1, adjCol1, r1, xMinusPA,
                    iIndex1, nA, jIndex, result);
        }
        else
        {
            /* New x' and sigma' is not feasible, so attempt
               a recovery. */
            getFeasiblePoint(getGamma2,  Gamma2,r2,
                             xMinusPB, iIndex2, violatedFace2,
                             result);

            /* Update (x - P ).
                            A
            */
            xMinusPA[0] = result[0] - p1[0];
            xMinusPA[1] = result[1] - p1[1];
            xMinusPA[2] = result[2] - p1[2];

            /* Update (x - P ).
                            B
            */
            xMinusPB[0] = result[0] - p2[0];
            xMinusPB[1] = result[1] - p2[1];
            xMinusPB[2] = result[2] - p2[2];

            /* Recover from infeasibility. */
            recover(getGamma1,  Gamma1,adjSurf1,
                    adjCol1, r1, xMinusPA, iIndex1,
                    nA, violatedFace1, 1, getGamma2, 
                    Gamma2, adjSurf2, adjCol2, r2, xMinusPB,
                    iIndex2, nB, jIndex, result);
        }
	  }
		/* Return appropriate code. */
		/*If the result of Sigma is negative, it means that the recovery fails. So, we need to
		return a restart signal.*/
		if(result[3] < TOLERANCE)
		{
			
			return 3;
		}
		else
		{
			return 1; 
		}
    
	}
    
}



/* ================================================================ */
/*                                                   T
   Function to compute the solution of the equation A x = b, where
    T
   A  is a 3 x 3 matrix.


   Input Parameters
   ================
   A[3][3] -- The 3 x 3 matrix "A".  Note, although the
              matrix "A" is passed in, this routine obtains the
                           T
              solution of A x = b, not Ax = b.
   b[3]    -- 3 x 1 vector containing the values for "b".


   Output Parameter
   ================
   x[3]    -- 3 x 1 vector to hold the result.


   Return Parameter
   ================
   A 0 indicates a solution is obtained.  Otherwise a 1 is returned
   to indicate that the matrix "A" given is singular.
*/
__inline int linearSolver3x3(double A[][3], double b[], double x[])
{
    /* Variables to store computed determinants. */
    double det012, det013, det032, det312;

    /* Array to hold computed cofactors so as to eliminate
       redundant computations. */
    double cfArray[4][4];

    /* Miscellaneous variables. */
    double overDet;


    /* ====================== */
    /* Determine the required cofactors. */
    cfArray[0][1] = A[1][0] * A[2][1] - A[2][0] * A[1][1];
    cfArray[0][2] = A[1][0] * A[2][2] - A[2][0] * A[1][2];
    cfArray[0][3] = A[1][0] * b[2] - A[2][0] * b[1];
    cfArray[1][2] = A[1][1] * A[2][2] - A[2][1] * A[1][2];
    cfArray[1][3] = A[1][1] * b[2] - A[2][1] * b[1];
    cfArray[3][2] = b[1] * A[2][2] - b[2] * A[1][2];

    det012 = A[0][0] * cfArray[1][2] - A[0][1] * cfArray[0][2] +
             A[0][2] * cfArray[0][1];

    if ((det012 > -TOLERANCE) && (det012 < TOLERANCE))
        /* Singular "A" matrix. */
        return 1;

    overDet = 1.0 / det012;

    det013 = A[0][0] * cfArray[1][3] - A[0][1] * cfArray[0][3] +
             b[0] * cfArray[0][1];
    det032 = A[0][0] * cfArray[3][2] - b[0] * cfArray[0][2] +
             A[0][2] * cfArray[0][3];
    det312 = b[0] * cfArray[1][2] - A[0][1] * cfArray[3][2] -
             A[0][2] * cfArray[1][3];

    x[0] = det312 * overDet;
    x[1] = det032 * overDet;
    x[2] = det013 * overDet;

    return 0;
}



/* ================================================================ */
/* In the case of edge-edge contact situation where feasibility is
   not maintained by either object, this function is called to obtain
   a new feasible point on one of the object. */
__inline void getFeasiblePoint(void (*getGamma)(), double *Gamma,double r[][3],
                               double xMinusP[], int iIndex[],
                               int violatedFace, double result[])
{
    /* Variables to hold temporary vectors. */
    double gammaA[3], gammaB[3], gammaC[3], v[3];

    /* Temporary variables. */
    double top, bottom, beta;


	/* ========================================== */
    /* Get direction vector formed by edge of active faces of
       infeasible object. */
    getGamma(Gamma, r, iIndex[0], gammaA);
    getGamma(Gamma, r, iIndex[1], gammaB);
    v[0] = gammaA[1] * gammaB[2] - gammaA[2] * gammaB[1];
    v[1] = gammaA[2] * gammaB[0] - gammaA[0] * gammaB[2];
    v[2] = gammaA[0] * gammaB[1] - gammaA[1] * gammaB[0];

    /* Compute "beta". */
    getGamma(Gamma, r, violatedFace, gammaC);
    top = result[3] - (gammaC[0] * xMinusP[0] +
                       gammaC[1] * xMinusP[1] +
                       gammaC[2] * xMinusP[2]);
    bottom = gammaC[0] * v[0] + gammaC[1] * v[1] + gammaC[2] * v[2];
    beta = top / bottom;

    /* Obtain new trial point. */
    result[0] += beta * v[0];
    result[1] += beta * v[1];
    result[2] += beta * v[2];
}



/* ================================================================ */
/* This function attempts to recover feasibility given that
   feasibility was not maintained. */
__inline void recover(void (*getGammaFailed)(), 
                      double *GammaFailed,int *adjSurfFailed,
                      int adjColFailed, double rFailed[][3],
                      double xMinusPFailed[], int iIndexFailed[],
                      int *nFailed, int violatedFace, int failedObject,
                      void (*getGamma)(), 
                      double *Gamma,int *adjSurf, int adjCol,
                      double r[][3], double xMinusP[], int iIndex[],
                      int *n, int jIndex[][2], double result[])
{
    /* Pointer to adjSurf arrays. */
    int *adjPtr;

    /* Variables to store temporary vectors. */
    double returnVal[3], deltaX[3], origin_deltaX[3], origin_xMinusP[3];
	double * Gamma_Ptr;

    /* Variable to hold active faces found. */
    int activeFace, activeFaceFailed;

    /* Miscellaneous variables. */
    double top, bottom, stepSize;
    double newSigma, value, reference;


    /* ========================================== */
    /* Need a reference value first. */
    activeFaceFailed = violatedFace;

	/*rotate the xMinusP to the original coordinate .*/
    get_origin_vector(xMinusPFailed, origin_xMinusP, rFailed); 
    
    /*                                       -1
          Compute "bottom" = Gamma .dx  . */
    Gamma_Ptr = GammaFailed + 3* activeFaceFailed;
	returnVal[0] = * Gamma_Ptr ++;
	returnVal[1] = * Gamma_Ptr ++;
	returnVal[2] = * Gamma_Ptr ;

    newSigma = returnVal[0] * origin_xMinusP[0] +
               returnVal[1] * origin_xMinusP[1] +
               returnVal[2] * origin_xMinusP[2];
    reference = newSigma ;
    
    /* Compute index into "adjSurfFailed". */
    adjPtr = adjSurfFailed + activeFaceFailed * adjColFailed;

    /* First we want to find the face belonging to the infeasible
       object that was pierced by the ray from the seed point of
       the infeasible object to the trial point. */
    do
    {
        /*                                       -1
              Compute "bottom" = Gamma .dx  . */
        Gamma_Ptr = GammaFailed + 3*(*adjPtr);
		returnVal[0] = * Gamma_Ptr ++;
		returnVal[1] = * Gamma_Ptr ++;
		returnVal[2] = * Gamma_Ptr ;

        value = returnVal[0] * origin_xMinusP[0] +
                returnVal[1] * origin_xMinusP[1] +
                returnVal[2] * origin_xMinusP[2];

        /* This algorithm employs the following strategy.  In
           searching through the adjacent faces for the
           greatest value, the first value greater than the
           current value will immediately result in branching
           to the associated surface. */
        if (value > reference)
        {
            newSigma = value;
            reference = value ;
            activeFaceFailed = *adjPtr;

            /* Move to new face. */
            adjPtr = adjSurfFailed + activeFaceFailed * adjColFailed;
        }
        else
            /* Next adjacent surface. */
            adjPtr++;
    } while (*adjPtr != -1);

    /* Get the new search direction. */
    deltaX[0] = -xMinusPFailed[0];
    deltaX[1] = -xMinusPFailed[1];
    deltaX[2] = -xMinusPFailed[2];

    /* Now let us do descend to obtain the second active 
       constraint from the feasible object. First, we need
       a reference value to kick off the search.  Begin with
       any face from active set of feasible object and compute
                            T
               "top" = gamma .deltaX - deltaSigma
                            i
                            T
                     = gamma .deltaX + newSigma
                            i
       since newSigma = -deltaSigma.
    */
    activeFace = iIndex[0];

	/*rotate the xMinusP to the original coordinate .*/
    get_origin_vector(deltaX, origin_deltaX ,r);
    get_origin_vector(xMinusP, origin_xMinusP,r); 
	
	/*                                       -1
          Compute "bottom" = Gamma .dx  . */
    Gamma_Ptr = Gamma + 3 * activeFace ;
	returnVal[0] = * Gamma_Ptr ++;
	returnVal[1] = * Gamma_Ptr ++;
	returnVal[2] = * Gamma_Ptr ;

 
    top = returnVal[0] * origin_deltaX[0] + returnVal[1] * origin_deltaX[1] +
          returnVal[2] * origin_deltaX[2] + newSigma;

    /*                                     T
       Compute "bottom" = (newSigma - gamma .(x - P))
    */
    bottom = newSigma - (returnVal[0] * origin_xMinusP[0] +
                         returnVal[1] * origin_xMinusP[1] +
                         returnVal[2] * origin_xMinusP[2]);

    /* Now obtain the reference step size. */
    stepSize = top / bottom;
    reference = stepSize ;

    /* Compute index into "adjSurf". */
    adjPtr = adjSurf + activeFace * adjCol;

    /* Here's our loop. */
    do
    {
		/*                                       -1
              Compute "bottom" = Gamma .dx  . */
        Gamma_Ptr = Gamma + 3 * (*adjPtr) ;
	    returnVal[0] = * Gamma_Ptr ++;
	    returnVal[1] = * Gamma_Ptr ++;
	    returnVal[2] = * Gamma_Ptr ;

 
        top = returnVal[0] * origin_deltaX[0] + 
			  returnVal[1] * origin_deltaX[1] +
              returnVal[2] * origin_deltaX[2] + newSigma;

        /*                                     T
           Compute "bottom" = (newSigma - gamma .(x - P))
        */
        bottom = newSigma - (returnVal[0] * origin_xMinusP[0] +
                             returnVal[1] * origin_xMinusP[1] +
                             returnVal[2] * origin_xMinusP[2]);


        /* Now obtain the step size. */
        value = top / bottom;

        /* This algorithm employs the following strategy.  In
           searching through the adjacent faces for the
           greatest value, the first value greater than the
           current value will immediately result in branching
           to the associated surface. */
        if (value > reference)
        {
            stepSize = value;
            reference = value ;
            activeFace = *adjPtr;

            /* Compute new index into "adjSurf". */
            adjPtr = adjSurf + activeFace * adjCol;
        }
        else
        {
            /* Next adjacent surface. */
            adjPtr++;
        }
    } while (*adjPtr != -1);

    /* Now that we have found the two active faces, we must update
       the indices. */
    *nFailed = *n = 1;
    iIndexFailed[0] = activeFaceFailed;
    iIndex[0] = activeFace;
    if (failedObject == 1)
    {
        jIndex[0][0] = activeFaceFailed;
        jIndex[0][1] = 1;
        jIndex[1][0] = activeFace;
        jIndex[1][1] = 2;
    }
    else
    {
        jIndex[0][0] = activeFace;
        jIndex[0][1] = 1;
        jIndex[1][0] = activeFaceFailed;
        jIndex[1][1] = 2;
    }

    /* Compute new "x" and "sigma" values. */
    stepSize = 1 / stepSize;
    result[0] += stepSize * deltaX[0];
    result[1] += stepSize * deltaX[1];
    result[2] += stepSize * deltaX[2];
    result[3] = (1.0 - stepSize) * newSigma;
}

