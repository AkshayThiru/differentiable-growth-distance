#include <stdio.h>
#include <limits.h>
#include <math.h>
#include "Gamma.h"
#include "Tolerance.h"
#include "incrCode4c.h"
#include "Solver4c.h"
#include "binarySearch.h"


/* Defines the threshold used to decide if the pseudo-binary search
   routine should be used. */
#define THRESHOLD   10


/* Local variables to store computed gamma values. */
double gammaVals[4][3];
double gammaAB, gammaAsq, gammaBsq, gammaABsq;


/* ================================================================ */
/* Local prototype declarations. */
__inline int initialize(double *Gamma1, double *Gamma2,
						int *adjSurf1, int *adjSurf2,
                        int adjCol1, int adjCol2, double r1[][3],
                        double r2[][3], double p1[], double p2[],int iIndex1[], int iIndex2[],
						int *nA, int *nB, int jIndex[][2], double result[], int increment);
__inline int faceFaceDescent(double *Gamma1, double *Gamma2,
                             int *adjSurf1, int *adjSurf2,
                             int *adjFaceCount1, int *adjFaceCount2,
                             int adjCol1, int adjCol2, double r1[][3],
                             double r2[][3], double p1[], double p2[],
                             int iIndex1[], int iIndex2[], int *nA,
                             int *nB, int jIndex[][2], double result[]);
__inline int faceEdgeDescent(double *Gamma1, double *Gamma2,
                             int *adjSurf1, int *adjSurf2,
                             int *adjFaceCount1, int *adjFaceCount2,
                             int adjCol1, int adjCol2, double r1[][3],
                             double r2[][3], double p1[], double p2[],
                             int iIndex1[], int iIndex2[], int *nA,
                             int *nB, int jIndex[][2], double result[]);
__inline int finalDescent(double *Gamma1, double *Gamma2,int *adjSurf1,
                          int *adjSurf2, int *adjFaceCount1,
                          int *adjFaceCount2, int adjCol1, int adjCol2,
                          double r1[][3], double r2[][3], double p1[],
                          double p2[], int iIndex1[], int iIndex2[],
                          int *nA, int *nB, int jIndex[][2],
                          double result[],double dInv[][4], int dInvReq);
int inIndex(int indexValue, int iIndexPtr[], int size);
__inline int linearSolver3x3c(double A[][3], int i0, int i1,
                              int i2, double x[]);




/* ================================================================ */
/* This is the entry point into this solver library.  It initializes
   all variables and prepares the necessary data required by the
   functions it calls.


   Input Parameters
   ================
  
   faces1,
   faces2          -- Integer values of the number of faces for objects 1 and 2
                      respectively.
   
   *Gamma1
   *Gamma2         -- Arrays of doubles containing the co-ordinates of the
                      adjusted outward normal vectors to the faces of objects 1 and
                      2 respectively. These values should be the output of the
                      routine setup_gamma. Note that setup_gamma is invoked only
					  once for each object. The main routine 'gd' for growth distance
					  computation calls to 'gd' use Gamma1 and Gamma2 directly. Refer to
                      documentations for setup_gamma for details.
  
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
   p2[3]         -- These 3 x 1 vectors specifices the translation
                    vectors associated with each object respectively.
                    Use these vectors to position the seed points
                    of the objects with respect to the global frame
                    of reference.
   iIndex1[3],
   iIndex2[3]    -- Index of active faces for objects 1 and 2.  These
                    two arrays are used internally by this routine.
                    This two indices should not be touched, in
                    particular, if incremental code is used.
   *nA, *nB      -- Count of number of elements in "iIndex1" and
                    "iIndex2" respectively.
   jIndex[4][2]  -- Contains indices NO. of active faces and which object
                    these active faces belongs to respectively.This 
					index is used internally.
   increment     -- Indicates if incremental codes should be used.

  
   OUTPUT PARAMETERS FROM GROWTH DISTANCE SOLVER
   =============================================
   result[4]       -- The first 3 x 1 vector (result[0] .. result[2])
                      contains the coordinates of the point of contact
                      of the objects corresponding to the sigma value
                      stored in "result[3]".
   dInv[4][4]      -- If "dInvReq" is set to 1, this 4 x 4 double array
                      returns the inverse of the C matrix.  See the appendix
					  below on how dInv can be used for gradient computations.



   Return Value
   ============
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
                  with 10000 faces.
*/
int solver(double *Gamma1, double * Gamma2,
		   int *adjSurf1, int *adjSurf2,
           int *adjFaceCount1, int *adjFaceCount2, int adjCol1,
           int adjCol2, double r1[][3], double r2[][3], double p1[],
           double p2[], int iIndex1[], int iIndex2[], int *nA, int *nB,
           int jIndex[][2],double result[], double dInv[][4], int dInvReq,int increment)
{

    int returnCode;


    if (increment == 1)
    {
        /* Use incremental code. */
        returnCode = incrSolver(Gamma1, Gamma2,
                                adjSurf1, adjSurf2, adjFaceCount1,
                                adjFaceCount2, adjCol1, adjCol2,
                                r1, r2, p1, p2, iIndex1, iIndex2, nA, nB,
                                jIndex, result);

#ifdef _CHECK_SEED
        /* Check for coincident seed points. */
        if (returnCode == 2)
		{ return returnCode;}
#endif

        /* Trap return codes. */
        if (returnCode == 0)
        {
            /* Feasibility has been maintained, so check for optimality. */
            returnCode = finalDescent(Gamma1, Gamma2,
									  adjSurf1, adjSurf2, adjFaceCount1,
                                      adjFaceCount2, adjCol1, adjCol2,
                                      r1, r2, p1, p2, iIndex1, iIndex2, nA,
                                      nB, jIndex, result, dInv, dInvReq);

            /*printf("Feasibility!");*/
			return returnCode;
        }

        if (returnCode == 1)
        {
            returnCode = faceFaceDescent(Gamma1, Gamma2,
				                         adjSurf1, adjSurf2, adjFaceCount1,
                                         adjFaceCount2, adjCol1, adjCol2,
                                         r1, r2, p1, p2, iIndex1, iIndex2, nA,
                                         nB, jIndex, result);

            
            if (returnCode != 0)
                return returnCode;

            returnCode = faceEdgeDescent(Gamma1,Gamma2,
										 adjSurf1, adjSurf2, adjFaceCount1,
                                         adjFaceCount2, adjCol1, adjCol2,
                                         r1, r2, p1, p2, iIndex1, iIndex2, nA,
                                         nB, jIndex, result);

            
            if (returnCode != 0)
                return returnCode;

            returnCode = finalDescent(Gamma1, Gamma2,
                                      adjSurf1, adjSurf2, adjFaceCount1,
                                      adjFaceCount2, adjCol1, adjCol2,
                                      r1, r2, p1, p2, iIndex1, iIndex2, nA,
                                      nB, jIndex, result, dInv, dInvReq);

            return returnCode;
        }

    }
    returnCode = initialize(Gamma1, Gamma2,
		                    adjSurf1, adjSurf2, adjCol1, adjCol2,
                            r1, r2, p1, p2, iIndex1, iIndex2, nA, nB,
                            jIndex,result,increment);

    /* Trap return codes. */
    if (returnCode != 0)
        return returnCode;

    returnCode = faceFaceDescent(Gamma1, Gamma2,
                                 adjSurf1, adjSurf2, adjFaceCount1,
                                 adjFaceCount2, adjCol1, adjCol2,
                                 r1, r2, p1, p2, iIndex1, iIndex2, nA, nB,
                                 jIndex, result);

    /* Trap return codes. */
    if (returnCode != 0)
        return returnCode;

    returnCode = faceEdgeDescent(Gamma1, Gamma2,
								 adjSurf1, adjSurf2, adjFaceCount1,
                                 adjFaceCount2, adjCol1, adjCol2,
                                 r1, r2, p1, p2, iIndex1, iIndex2, nA, nB,
                                 jIndex, result);

    /* Trap return codes. */
    if (returnCode != 0)
        return returnCode;

    returnCode = finalDescent(Gamma1, Gamma2,
							  adjSurf1, adjSurf2, adjFaceCount1,
                              adjFaceCount2, adjCol1, adjCol2,
                              r1, r2, p1, p2, iIndex1, iIndex2, nA, nB,
                              jIndex, result,dInv, dInvReq);

    return returnCode;
}



/* ================================================================ */
/* Function to find the two faces punctured by the direction
   vector joining the two seed points.  With this, we will be able
   to obtain the first two active constraints.


   Input Parameters
   ================
   As in function "solver".


   Input/Output Parameter
   ======================
   result[4] -- The first 3 x 1 vector (result[0] .. result[2])
                contains the x coordinate of the point of contact
                of the objects corresponding to the sigma value
                stored in "result[3]".


   Return Value
   ============
   The function returns a 1 if the values contained in the "result"
   variable is the solution of the problem, otherwise a 0 is
   returned.
*/
__inline int initialize(double *Gamma1,double *Gamma2,
						int *adjSurf1, int *adjSurf2,
                        int adjCol1, int adjCol2, double r1[][3],
                        double r2[][3], double p1[], double p2[],int iIndex1[], int iIndex2[],
						int *nA, int *nB,int jIndex[][2], double result[], int increment)
{
    /* Pointer to adjSurf arrays. */
    int *adjPtr;

	/*array store the negative of P2-P1 to make calling begin)face funnction convinent*/
	double dp[3],  minusdp[3];
    
	/* Miscellaneous variables. */
    int surface1, surface2;
    double alpha1, alpha2, reference;
    double lamda, value, value1, value2;


    /* ========================================== */
#ifdef CHECK_SEED
    /* Get direction vector joining the two seed points. */
    dP[0] = p2[0] - p1[0];
    dP[1] = p2[1] - p1[1];
    dP[2] = p2[2] - p1[2];


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

	dp[0] = p2[0] - p1[0];
    dp[1] = p2[1] - p1[1];
    dp[2] = p2[2] - p1[2];
	minusdp[0]=p1[0]-p2[0];
	minusdp[1]=p1[1]-p2[1];
	minusdp[2]=p1[2]-p2[2];
    
    /*set initial faces.*/
	surface1 = 0;   
	surface2 = 0;   

    /*           T
                n    T
       Compute ----.R .(p2 - p1).
                T -
               n .x
    */
   
	value1 = getGamma1P1(Gamma1, r1, p1, surface1);
    value2 = getGamma1P2(Gamma1, r1, p2, surface1);
    alpha1 = value2 - value1;
    reference = alpha1 ;


    /* Compute initial index into "adjSurf1". */
    adjPtr = adjSurf1 + surface1 * adjCol1;


    /* Check only adjacent surfaces. */
    do
    {
        /*           T
                    n    T
           Compute ----.R .(p2 - p1).
                    T -
                   n .x
        */
        value1 = getGamma1P1(Gamma1, r1, p1, *adjPtr);
        value2 = getGamma1P2(Gamma1, r1, p2, *adjPtr);
        value = value2 - value1;


        /* This algorithm employs the following strategy.  In
           searching through the adjacent surfaces for the
           greatest gamma value, the first value greater value
           encountered will immediately result in branching
           to the associated surface. */
        if (value > reference)
        {
            alpha1 = value;
            surface1 = *adjPtr;
            reference = alpha1 ;


            /* Compute new index into "adjSurf1". */
            adjPtr = adjSurf1 + surface1 * adjCol1;

        }
        else
            /* Next adjacent surface. */
            adjPtr++;
    } while (*adjPtr != -1);


    /* Update "iIndex1". */
    iIndex1[0] = surface1;
    *nA = 1;

    /* Update jIndex. */
    jIndex[0][0] = surface1;
    jIndex[0][1] = 1;


    /* ========================================== */
    /*           T
                n    T
       Compute ----.R .(p1 - p2).
                T -
               n .x
    */
	value1 = getGamma2P1(Gamma2, r2, p1, surface2);
    value2 = getGamma2P2(Gamma2, r2, p2, surface2);
    alpha2 = value1 - value2;
    reference = alpha2 ;


    /* Compute initial index into "adjSurf2". */
    adjPtr = adjSurf2 + surface2 * adjCol2;


    /* Check only adjacent surfaces. */
    do
    {
        /*           T
                    n    T
           Compute ----.R .(p1 - p2).
                    T -
                   n .x
        */
        value1 = getGamma2P1(Gamma2, r2, p1, *adjPtr);
        value2 = getGamma2P2(Gamma2, r2, p2, *adjPtr);
        value = value1 - value2;

		  /* This algorithm employs the following strategy.  In
           searching through the adjacent surfaces for the
           greatest gamma value, the first value greater value
           encountered will immediately result in branching
           to the associated surface. */
        if (value > reference)
        {
            alpha2 = value;
            surface2 = *adjPtr;
            reference = alpha2 ;



            /* Compute new index into "adjSurf2". */
            adjPtr = adjSurf2 + surface2 * adjCol2;

        }
        else
            /* Next adjacent surface. */
            adjPtr++;
    } while (*adjPtr != -1);


    /* Update "iIndex2". */
    iIndex2[0] = surface2;
    *nB = 1;

    /* Update "jIndex". */
    jIndex[1][0] = surface2;
    jIndex[1][1] = 2;

    /* Compute result of run. */
    result[3] = (alpha1 * alpha2) / (alpha1 + alpha2);
    lamda = 1.0 / alpha2 * result[3];
    result[0] = p2[0] + lamda * (p1[0] - p2[0]);
    result[1] = p2[1] + lamda * (p1[1] - p2[1]);
    result[2] = p2[2] + lamda * (p1[2] - p2[2]);

    return 0;
}



/* ================================================================ */
/* Function to do a face-face descent given 2 active constraints.


   Input Parameters
   ================
   As in function "solver".


   Input/Output Parameter
   ======================
   result[4] -- The first 3 x 1 vector (result[0] .. result[2])
                contains the x coordinate of the point of contact
                of the objects corresponding to the sigma value
                stored in "result[3]".


   Return Value
   ============
   The function returns a 1 if the values contained in the "result"
   variable is the solution of the problem, otherwise a 0 is
   returned.
*/
__inline int faceFaceDescent(double *Gamma1, double *Gamma2,
                             int *adjSurf1, int *adjSurf2,
                             int *adjFaceCount1, int *adjFaceCount2,
                             int adjCol1, int adjCol2, double r1[][3],
                             double r2[][3], double p1[], double p2[],
                             int iIndex1[], int iIndex2[], int *nA,
                             int *nB, int jIndex[][2], double result[])
{
    /* Pointer to adjSurf arrays. */
    int *adjPtr;

    /* Miscellaneous variables. */
    int activeFace1, activeFace2;
    double stepSize1, stepSize2;
    double top, bottom, overBottom, reference, value, testCondition;
    double betaA, betaB, returnVal[3];
	double deltaX[3], origin_deltaX[3], *Gamma_Ptr ;
    double xMinusP[3], origin_xMinusP[3];
 
    /* ========================================== */
    /* Compute "bottom". */
    getGamma1(Gamma1, r1, iIndex1[0], gammaVals[0]);
    getGamma2(Gamma2, r2, iIndex2[0], gammaVals[1]);

    gammaAsq = gammaVals[0][0] * gammaVals[0][0] +
               gammaVals[0][1] * gammaVals[0][1] +
               gammaVals[0][2] * gammaVals[0][2];
    gammaBsq = gammaVals[1][0] * gammaVals[1][0] +
               gammaVals[1][1] * gammaVals[1][1] +
               gammaVals[1][2] * gammaVals[1][2];
    gammaAB = gammaVals[0][0] * gammaVals[1][0] +
              gammaVals[0][1] * gammaVals[1][1] +
              gammaVals[0][2] * gammaVals[1][2];
    gammaABsq = gammaAB * gammaAB;
	/*testCondition = 1- cos( face(a) - face(b)). If the value of the cos is 1, it means 
	the two faces are parallel.*/
	testCondition = 1 - gammaABsq/(gammaAsq*gammaBsq);
    bottom = gammaAsq * gammaBsq - gammaABsq;


    /* See if optimality has been reached. */
    if (testCondition < TOLERANCE)
	{ jIndex[2][0]=-1;jIndex[2][1]=-1;jIndex[3][0]=-1;jIndex[3][1]=-1; return 1;}

    /* Invert bottom so that we can use multiplication
       instead of division (which is less efficient). */
    overBottom = 1.0 / bottom;

    /* Compute "betaA". */
    betaA = (gammaAB - gammaBsq) * overBottom;

    /* Compute "betaB". */
    betaB = (gammaAB - gammaAsq) * overBottom;

    /* Compute "deltaX". */
    deltaX[0] = betaA * gammaVals[0][0] + betaB * gammaVals[1][0];
    deltaX[1] = betaA * gammaVals[0][1] + betaB * gammaVals[1][1];
    deltaX[2] = betaA * gammaVals[0][2] + betaB * gammaVals[1][2];

    /* Store (x    - P ).
               init   A
    */
    xMinusP[0] = result[0] - p1[0];
    xMinusP[1] = result[1] - p1[1];
    xMinusP[2] = result[2] - p1[2];

    /* Compute index into "adjSurf1". */
    adjPtr = adjSurf1 + iIndex1[0] * adjCol1;

    /* Find largest feasible step size. */
    if (adjFaceCount1[iIndex1[0]] > THRESHOLD)
    {
        binarySearch(getGamma1,Gamma1, r1, deltaX, xMinusP,
                     result[3], adjPtr, adjFaceCount1[iIndex1[0]],
                     &activeFace1, &stepSize1);
        reference = stepSize1 ;

        /* Set "setSize2" to a ridiculously large value. */
        stepSize2 = 1e10;
    }
    else
    {
        /* Set "setSize.." to a ridiculously large value. */
        reference = stepSize1 = stepSize2 = 1e10;

        /*rotate the deltaX and xMinusP to the original coordinate .*/
        get_origin_vector(deltaX, origin_deltaX ,r1);
        get_origin_vector(xMinusP, origin_xMinusP, r1); 
        do
        {
            /*                         -1
               Compute "bottom" = Gamma .dx  . */
            Gamma_Ptr = Gamma1 + 3*(*adjPtr);
			returnVal[0] = * Gamma_Ptr ++;
			returnVal[1] = * Gamma_Ptr ++;
			returnVal[2] = * Gamma_Ptr ;
            			
            bottom = returnVal[0] * origin_deltaX[0] +
                     returnVal[1] * origin_deltaX[1] +
                     returnVal[2] * origin_deltaX[2];


            if (bottom >  - 1.0)
            {
                /*                           T
                   Compute (sigma     - gamma  .(x     - P )).
                                 init        ic   init    A
                */
                top = result[3] - (returnVal[0] * origin_xMinusP[0] +
                                   returnVal[1] * origin_xMinusP[1] +
                                   returnVal[2] * origin_xMinusP[2]);


                /* Computed value of "alpha3" must be positive,
                   otherwise the value of growth factor "sigma"
                   will increase. */
                if (top > 0.0)
                {
                    value = top / (bottom + 1.0);


                    if (value < reference)
                    {
                        stepSize1 = value;
                        activeFace1 = *adjPtr;
                        reference = value ;

                    }
                }
            }

            /* Next adjacent surface. */
            adjPtr++;
        } while (*adjPtr != -1);
    }

    /* ========================================== */
    /* Store (x    - P ).
               init   B
    */
    xMinusP[0] = result[0] - p2[0];
    xMinusP[1] = result[1] - p2[1];
    xMinusP[2] = result[2] - p2[2];

    /* Compute index into "adjSurf2". */
    adjPtr = adjSurf2 + iIndex2[0] * adjCol2;

    /* Find largest feasible step size. */
    if (adjFaceCount2[iIndex2[0]] > THRESHOLD)
        binarySearch(getGamma2, Gamma2, r2, deltaX, xMinusP,
                     result[3], adjPtr, adjFaceCount2[iIndex2[0]], &activeFace2,
                     &stepSize2);
    else
    {
        
		/*rotate the deltaX and xMinusP to the original coordinate .*/
        get_origin_vector(deltaX, origin_deltaX ,r2);
        get_origin_vector(xMinusP, origin_xMinusP, r2); 
		do
        {
            
		 /*                                    -1
               Compute "bottom" = Gamma .dx  . */
            Gamma_Ptr = Gamma2 + 3*(*adjPtr);
			returnVal[0] = * Gamma_Ptr ++;
			returnVal[1] = * Gamma_Ptr ++;
			returnVal[2] = * Gamma_Ptr ;
	
            bottom = returnVal[0] * origin_deltaX[0] +
                     returnVal[1] * origin_deltaX[1] +
                     returnVal[2] * origin_deltaX[2];

            if (bottom >  - 1.0)
            {
                /*                           T
                   Compute (sigma     - gamma  .(x     - P )).
                                 init        ic   init    B
                */
                top = result[3] - (returnVal[0] * origin_xMinusP[0] +
                                   returnVal[1] * origin_xMinusP[1] +
                                   returnVal[2] * origin_xMinusP[2]);


                /* Computed value of "alpha3" must be positive,
                   otherwise the value of growth factor "sigma"
                   will increase. */
                if (top > 0.0)
                {
                    value = top / (bottom + 1.0);


                    if (value < reference)
                    {
                        stepSize2 = value;
                        activeFace2 = *adjPtr;
                        reference = value ;

                    }
                }
            }

            /* Next adjacent surface. */
            adjPtr++;
        } while (*adjPtr != -1);
    }

    /* Decide which is the active face is appropriate. */
    if (stepSize1 < stepSize2)
    {
        /* Update "iIndex1". */
        iIndex1[1] = activeFace1;
        (*nA)++;

        /* Update "jIndex". */
        jIndex[2][0] = activeFace1;
        jIndex[2][1] = 1;

        /* Compute result of run. */
        result[0] += stepSize1 * deltaX[0];
        result[1] += stepSize1 * deltaX[1];
        result[2] += stepSize1 * deltaX[2];
        result[3] -= stepSize1;
    }
    else
    {
        /* Update "iIndex2". */
        iIndex2[1] = activeFace2;
        (*nB)++;

        /* Update "jIndex". */
        jIndex[2][0] = activeFace2;
        jIndex[2][1] = 2;

        /* Compute result of run. */
        result[0] += stepSize2 * deltaX[0];
        result[1] += stepSize2 * deltaX[1];
        result[2] += stepSize2 * deltaX[2];
        result[3] -= stepSize2;
    }

    return 0;
}

/* ================================================================ */
/* Function to do a face-edge descent given 3 active constraints.


   Input Parameters
   ================
   As in function "solver".


   Input/Output Parameter
   ======================
   result[4] -- The first 3 x 1 vector (result[0] .. result[2])
                contains the x coordinate of the point of contact
                of the objects corresponding to the sigma value
                stored in "result[3]".


   Return Value
   ============
   The return values of this function may be interpreted as follows:
        0:    No errors.
        1:    Solution is optimal.
        2:    Divide by zero error.  Possible error in definitions of
              input data.
*/
__inline int faceEdgeDescent(double *Gamma1, double *Gamma2,
                             int *adjSurf1, int *adjSurf2,
                             int *adjFaceCount1, int *adjFaceCount2,
                             int adjCol1, int adjCol2, double r1[][3],
                             double r2[][3], double p1[], double p2[],
                             int iIndex1[], int iIndex2[], int *nA,
                             int *nB, int jIndex[][2], double result[])
{
    /* Pointer to adjSurf arrays. */
    int *adjPtr;

    /* Variable to indicate which object the third constraint
       was obtained from. */
    int whichIndex;

    /* Miscellaneous variables. */
    int surface4;
    double alpha4, value, top, bottom, reference;
    double returnVal[3];
	double deltaX[3], origin_deltaX[3];
	double xMinusP[3], origin_xMinusP[3];
	double  * Gamma_Ptr;


    /* ========================================== */
    /* Which object did the third constraint come from? */
    if (jIndex[2][1] == 1)
    {
        /* Get gamma values for this new constraint first. */
        getGamma1(Gamma1, r1, iIndex1[1], gammaVals[2]);

        /* Compute search direction "deltaX". */
        if (linearSolver3x3c(gammaVals, 0, 1, 2, deltaX) == 1)
            /* Singular matrix.  Solution is optimal. */
		{jIndex[3][0]=-1;jIndex[3][1]=-1;return 1;}

        /* Store (x    - P ).
                   init   B
        */
        xMinusP[0] = result[0] - p2[0];
        xMinusP[1] = result[1] - p2[1];
        xMinusP[2] = result[2] - p2[2];

        /* Check adjacent surfaces of active face from object 2. */
        adjPtr = adjSurf2 + iIndex2[0] * adjCol2;

        if (adjFaceCount2[iIndex2[0]] > THRESHOLD)
        {
            binarySearch(getGamma2, Gamma2, r2, deltaX, xMinusP,
                         result[3], adjPtr, adjFaceCount2[iIndex2[0]],
                         &surface4, &alpha4);
            reference = alpha4 ;
            whichIndex = 2;
        }
        else
        {
            /* Set "alpha4" to a ridiculously large value. */
            reference = alpha4 = 1e10;

            /*rotate the deltaX and xMinusP to the original coordinate .*/
            get_origin_vector(deltaX, origin_deltaX ,r2);
            get_origin_vector(xMinusP, origin_xMinusP, r2);       
			
			
			do
            {
			    /*                                 -1
                Compute "bottom" = Gamma .dx  . */
                Gamma_Ptr = Gamma2 + 3*(*adjPtr);
			    returnVal[0] = * Gamma_Ptr ++;
			    returnVal[1] = * Gamma_Ptr ++;
			    returnVal[2] = * Gamma_Ptr ;

                bottom = returnVal[0] * origin_deltaX[0] +
                         returnVal[1] * origin_deltaX[1] +
                         returnVal[2] * origin_deltaX[2];


                if (bottom >  - 1.0)
                {
                    /*                           T
                       Compute (sigma     - gamma  .(x    - P )).
                                     init        ic   init   B
                    */
                    top = result[3] - (returnVal[0] * origin_xMinusP[0] +
                                       returnVal[1] * origin_xMinusP[1] +
                                       returnVal[2] * origin_xMinusP[2]);


                    /* Computed value of "alpha4" must be positive,
                       otherwise the value of growth factor "sigma"
                       will increase. */
                    if (top > 0.0)
                    {
                        value = top / (bottom + 1.0);


                        if (value < reference)
                        {
                            alpha4 = value;
                            surface4 = *adjPtr;
                            reference = alpha4 ;

                            /* Make sure we know where we got our "alpha4"
                               value. */
                            whichIndex = 2;
                        }
                    }
                }

                /* Next adjacent surface. */
                adjPtr++;
            } while (*adjPtr != -1);
        }

        /* Store (x    - P ).
                   init   A
        */
        xMinusP[0] = result[0] - p1[0];
        xMinusP[1] = result[1] - p1[1];
        xMinusP[2] = result[2] - p1[2];

		/*rotate the deltaX and xMinusP to the original coordinate .*/
        get_origin_vector(deltaX, origin_deltaX ,r1);
        get_origin_vector(xMinusP, origin_xMinusP, r1); 

        /* Third constraint is from object 1.  We need to check
           adjacent faces of the active face from object 1 that
           has the lower number of adjacent faces. */
        if ((adjFaceCount1[iIndex1[0]]) < (adjFaceCount1[iIndex1[1]]))
        {
            adjPtr = adjSurf1 + iIndex1[0] * adjCol1;

            do
            {
                /* Must make sure the active faces are not included. */
                if (*adjPtr != iIndex1[1])
                {
                    /*                                    -1
                       Compute "bottom" = Gamma .dx  . */
                    Gamma_Ptr = Gamma1 + 3*(*adjPtr);
			        returnVal[0] = * Gamma_Ptr ++;
			        returnVal[1] = * Gamma_Ptr ++;
			        returnVal[2] = * Gamma_Ptr ;
                    
				    bottom = returnVal[0] * origin_deltaX[0] +
                             returnVal[1] * origin_deltaX[1] +
                             returnVal[2] * origin_deltaX[2];

                    if (bottom >  - 1.0)
                    {
                        /*                           T
                           Compute (sigma     - gamma  .(x    - P )).
                                         init        ic   init   A
                        */
                        top = result[3] - (returnVal[0] * origin_xMinusP[0] +
                                           returnVal[1] * origin_xMinusP[1] +
                                           returnVal[2] * origin_xMinusP[2]);


                        /* Computed value of "alpha4" must be positive,
                           otherwise the value of growth factor "sigma"
                           will increase. */
                        if (top > 0.0)
                        {
                            value = top / (bottom + 1.0);


                            if (value < reference)
                            {
                                alpha4 = value;
                                surface4 = *adjPtr;
                                reference = alpha4 ;

                                /* Make sure we know where we got our "alpha4"
                                   value. */
                                whichIndex = 1;
                            }
                        }
                    }
                }

                /* Next adjacent surface. */
                adjPtr++;
            } while (*adjPtr != -1);
        }
        else
        {
            adjPtr = adjSurf1 + iIndex1[1] * adjCol1;

            do
            {
                /* Must make sure the active faces are not included. */
                if (*adjPtr != iIndex1[0])
                {
                    /*                                    -1
                       Compute "bottom" = Gamma .dx  . */
                    Gamma_Ptr = Gamma1 + 3*(*adjPtr);
			        returnVal[0] = * Gamma_Ptr ++;
			        returnVal[1] = * Gamma_Ptr ++;
			        returnVal[2] = * Gamma_Ptr ;
                    
				    bottom = returnVal[0] * origin_deltaX[0] +
                             returnVal[1] * origin_deltaX[1] +
                             returnVal[2] * origin_deltaX[2];

                    if (bottom >  - 1.0)
                    {
                        /*                           T
                           Compute (sigma     - gamma  .(x    - P )).
                                         init        ic   init   A
                        */
                        top = result[3] - (returnVal[0] * origin_xMinusP[0] +
                                           returnVal[1] * origin_xMinusP[1] +
                                           returnVal[2] * origin_xMinusP[2]);

                        /* Computed value of "alpha4" must be positive,
                           otherwise the value of growth factor "sigma"
                           will increase. */
                        if (top > 0.0)
                        {
                            value = top / (bottom + 1.0);

                            if (value < reference)
                            {
                                alpha4 = value;
                                surface4 = *adjPtr;
                                reference = alpha4 ;


                                /* Make sure we know where we got our "alpha4"
                                   value. */
                                whichIndex = 1;
                            }
                        }
                    }
                }

                /* Next adjacent surface. */
                adjPtr++;
            } while (*adjPtr != -1);
        }
    }
    else
    {
        /* Get gamma values for this new constraint first. */
        getGamma2(Gamma2, r2, iIndex2[1], gammaVals[2]);

        /* Compute search direction "deltaX". */
        if (linearSolver3x3c(gammaVals, 0, 1, 2, deltaX) == 1)
            /* Singular matrix.  Solution is optimal. */
		{  jIndex[3][0]=-1;jIndex[3][1]=-1; return 1;}

        /* Store (x    - P ).
                   init   A
        */
        xMinusP[0] = result[0] - p1[0];
        xMinusP[1] = result[1] - p1[1];
        xMinusP[2] = result[2] - p1[2];

        /* Third constraint is from object 2.  We start checking
           adjacent surfaces of active face from object 1. */
        adjPtr = adjSurf1 + iIndex1[0] * adjCol1;

        if (adjFaceCount1[iIndex1[0]] > THRESHOLD)
        {
            binarySearch(getGamma1, Gamma1, r1, deltaX, xMinusP,
                         result[3], adjPtr, adjFaceCount1[iIndex1[0]],
                         &surface4, &alpha4);
            reference = alpha4 ;
            whichIndex = 1;
        }
        else
        {
            /* Set "alpha4" to a ridiculously large value. */
            reference = alpha4 = 1e10;

	        /*rotate the deltaX and xMinusP to the original coordinate .*/
            get_origin_vector(deltaX, origin_deltaX ,r1);
            get_origin_vector(xMinusP, origin_xMinusP, r1); 

            do
            {
                /*                                    -1
                   Compute "bottom" = Gamma .dx  . */
                Gamma_Ptr = Gamma1 + 3*(*adjPtr);
			    returnVal[0] = * Gamma_Ptr ++;
			    returnVal[1] = * Gamma_Ptr ++;
			    returnVal[2] = * Gamma_Ptr ;
	
                bottom = returnVal[0] * origin_deltaX[0] +
                         returnVal[1] * origin_deltaX[1] +
                         returnVal[2] * origin_deltaX[2];

                if (bottom > - 1.0)
                {
                    /*                           T
                       Compute (sigma     - gamma  .(x    - P )).
                                     init        ic   init   A
                    */
                    top = result[3] - (returnVal[0] * origin_xMinusP[0] +
                                       returnVal[1] * origin_xMinusP[1] +
                                       returnVal[2] * origin_xMinusP[2]);

                    /* Computed value of "alpha4" must be positive,
                       otherwise the value of growth factor "sigma"
                       will increase. */
                    if (top > 0)
                    {
                        value = top / (bottom + 1.0);

                        if (value < reference)
                        {
                            alpha4 = value;
                            surface4 = *adjPtr;
                            reference = alpha4 ;

                            /* Make sure we know where we got our "alpha4"
                               value. */
                            whichIndex = 1;
                        }
                    }
                }

                /* Next adjacent surface. */
                adjPtr++;
            } while (*adjPtr != -1);
        }

        /* Store (x    - P ).
                   init   B
        */
        xMinusP[0] = result[0] - p2[0];
        xMinusP[1] = result[1] - p2[1];
        xMinusP[2] = result[2] - p2[2];

        /*rotate the deltaX and xMinusP to the original coordinate .*/
        get_origin_vector(deltaX, origin_deltaX ,r2);
        get_origin_vector(xMinusP, origin_xMinusP, r2); 

		/* Check adjacent surfaces of active face from object 2 that
           has the lower number of adjacent faces. */
        if (adjFaceCount2[iIndex2[0]] < adjFaceCount2[iIndex2[1]])
        {
            adjPtr = adjSurf2 + iIndex2[0] * adjCol2;

            do
            {
                /* Must make sure the active faces are not included. */
                if (*adjPtr != iIndex2[1])
                {
                    /*                                    -1
                       Compute "bottom" = Gamma .dx  . */
                    Gamma_Ptr = Gamma2 + 3*(*adjPtr);
			        returnVal[0] = * Gamma_Ptr ++;
			        returnVal[1] = * Gamma_Ptr ++;
			        returnVal[2] = * Gamma_Ptr ;

                    bottom = returnVal[0] * origin_deltaX[0] +
                             returnVal[1] * origin_deltaX[1] +
                             returnVal[2] * origin_deltaX[2];

                    if (bottom >  - 1.0)
                    {
                        /*                           T
                           Compute (sigma     - gamma  .(x    - P )).
                                         init        ic   init   B
                        */
                        top = result[3] - (returnVal[0] * origin_xMinusP[0] +
                                           returnVal[1] * origin_xMinusP[1] +
                                           returnVal[2] * origin_xMinusP[2]);

                        /* Computed value of "alpha4" must be positive,
                           otherwise the value of growth factor "sigma"
                           will increase. */
                        if (top > 0.0)
                        {
                            value = top / (bottom + 1.0);

                            if (value < reference)
                            {
                                alpha4 = value;
                                surface4 = *adjPtr;
                                reference = alpha4 ;

                                /* Make sure we know where we got our "alpha4"
                                   value. */
                                whichIndex = 2;
                            }
                        }
                    }
                }

                /* Next adjacent surface. */
                adjPtr++;
            } while (*adjPtr != -1);
        }
        else
        {
            adjPtr = adjSurf2 + iIndex2[1] * adjCol2;

            do
            {
                /* Must make sure the active faces are not included. */
                if (*adjPtr != iIndex2[0])
                {
                    /*                                    -1
                       Compute "bottom" = Gamma .dx  . */
                    Gamma_Ptr = Gamma2 + 3*(*adjPtr);
			        returnVal[0] = * Gamma_Ptr ++;
			        returnVal[1] = * Gamma_Ptr ++;
			        returnVal[2] = * Gamma_Ptr ;
                    
                    bottom = returnVal[0] * origin_deltaX[0] +
                             returnVal[1] * origin_deltaX[1] +
                             returnVal[2] * origin_deltaX[2];

                    if (bottom > - 1.0)
                    {
                        /*                           T
                           Compute (sigma     - gamma  .(x    - P )).
                                         init        ic   init   B
                        */
                        top = result[3] - (returnVal[0] * origin_xMinusP[0] +
                                           returnVal[1] * origin_xMinusP[1] +
                                           returnVal[2] * origin_xMinusP[2]);

                        /* Computed value of "alpha4" must be positive,
                           otherwise the value of growth factor "sigma"
                           will increase. */
                        if (top > 0.0)
                        {
                            value = top / (bottom + 1.0);

                            if (value < reference)
                            {
                                alpha4 = value;
                                surface4 = *adjPtr;
                                reference = alpha4 ;

                                /* Make sure we know where we got our "alpha4"
                                   value. */
                                whichIndex = 2;
                            }
                        }
                    }
                }

                /* Next adjacent surface. */
                adjPtr++;
            } while (*adjPtr != -1);
        }
    }

    /* Obtain "gGamma" values for new constraint. */
    if (whichIndex == 1)
    {
        getGamma1(Gamma1, r1, surface4, gammaVals[3]);

        /* Update "iIndex1". */
        iIndex1[*nA] = surface4;
        (*nA)++;
    }
    else
    {
        getGamma2(Gamma2, r2, surface4, gammaVals[3]);

        /* Update "iIndex2". */
        iIndex2[*nB] = surface4;
        (*nB)++;
    }

    /* Update "jIndex". */
    jIndex[3][0] = surface4;
    jIndex[3][1] = whichIndex;

    /* Compute result of run. */
    result[0] += alpha4 * deltaX[0];
    result[1] += alpha4 * deltaX[1];
    result[2] += alpha4 * deltaX[2];
    result[3] -= alpha4;

    return 0;
}

/* ================================================================ */
/* Function to do a edge-edge or a face-vertex descent given 4
   active constraints.


   Input Parameters
   ================
   As in function "solver".


   Input/Output Parameter
   ======================
   result[4] -- The first 3 x 1 vector (result[0] .. result[2])
                contains the x coordinate of the point of contact
                of the objects corresponding to the sigma value
                stored in "result[3]".


   Return Value
   ============
   The return values of this function may be interpreted as follows:
        1:    Solution has been reached.  The values contained in
              the "result" variable is the solution to the problem.
        2:    Divide by zero error.  Possible error in definitions of
              input data.
*/
__inline int finalDescent(double *Gamma1, double *Gamma2, 
						  int *adjSurf1,int *adjSurf2, int *adjFaceCount1,
                          int *adjFaceCount2, int adjCol1, int adjCol2,
                          double r1[][3], double r2[][3], double p1[],
                          double p2[], int iIndex1[], int iIndex2[],
                          int *nA, int *nB, int jIndex[][2],
                          double result[],double dInv[][4], int dInvReq)
{
    /* Pointer to adjSurf arrays. */
    int *adjPtr;

    /* Variable to indicate from which object the third constraint
       was obtained from. */
    int whichIndex;

    /* Variable to find search direction. */
    int dIndex;
    double descent;

    /* Laranage multipliers. */
    double lamda[4];

  /* Array to store cofactors used in determining determinants. */
    double cfArray[4][4];
    double cf1, cf2, cf3;
	double cf1000, cf2000, cf2010, cf3000, cf3010, cf3020;

    /* Variables to hold the result of the determinant computations. */
    double det012, det013, det032, det312;
    double overDet;

    /* Variable to check for optimality. */
    int optimal;

    /* Miscellaneous variables. */
    register int i, k;
    int surface, obj1Face, obj2Face;
    int obj1AdjFaceCount, obj2AdjFaceCount;
    int obj1ActiveFaceCount, obj2ActiveFaceCount;
    double alpha, value, top, bottom, reference, Max_det, testCondition;
    double returnVal[3], deltaX[3], origin_deltaX[3];
	double xMinusPA[3], origin_xMinusPA[3];
	double xMinusPB[3], origin_xMinusPB[3];
	double * Gamma_Ptr;


    /* ========================================== */
    /* Compute Lagrange multipliers.  */
    cfArray[0][1] = gammaVals[0][1] * gammaVals[1][2] -
                    gammaVals[0][2] * gammaVals[1][1];
    cfArray[0][2] = gammaVals[0][1] * gammaVals[2][2] -
                    gammaVals[0][2] * gammaVals[2][1];
    cfArray[0][3] = gammaVals[0][1] * gammaVals[3][2] -
                    gammaVals[0][2] * gammaVals[3][1];
    cfArray[1][2] = gammaVals[1][1] * gammaVals[2][2] -
                    gammaVals[1][2] * gammaVals[2][1];
    cfArray[1][3] = gammaVals[1][1] * gammaVals[3][2] -
                    gammaVals[1][2] * gammaVals[3][1];
    cfArray[3][2] = gammaVals[3][1] * gammaVals[2][2] -
                    gammaVals[3][2] * gammaVals[2][1];

    det012 = gammaVals[0][0] * cfArray[1][2] -
             gammaVals[1][0] * cfArray[0][2] +
             gammaVals[2][0] * cfArray[0][1];
    det013 = gammaVals[0][0] * cfArray[1][3] -
             gammaVals[1][0] * cfArray[0][3] +
             gammaVals[3][0] * cfArray[0][1];
    det032 = gammaVals[0][0] * cfArray[3][2] -
             gammaVals[3][0] * cfArray[0][2] +
             gammaVals[2][0] * cfArray[0][3];
    det312 = gammaVals[3][0] * cfArray[1][2] -
             gammaVals[1][0] * cfArray[3][2] -
             gammaVals[2][0] * cfArray[1][3];
    bottom = det012 - det013 - det032 - det312;

    /*calculate test condition.*/
	Max_det = fabs(det012);
	if ( fabs(det013) > Max_det) Max_det = fabs(det013);
	if ( fabs(det032) > Max_det) Max_det = fabs(det032);	
	if ( fabs(det312) > Max_det) Max_det = fabs(det312);
    testCondition = fabs(bottom)/Max_det;

    /* Check for divide-by-zero. */
    if (testCondition < TOLERANCE)
        return 3;

    overDet = 1.0 / bottom;

    lamda[0] = -(det312 * overDet);
    lamda[1] = -(det032 * overDet);
    lamda[2] = -(det013 * overDet);
    lamda[3] = 1.0 - (lamda[0] + lamda[1] + lamda[2]);

    for (;;)
    {
        /* Assume optimal. */
        optimal = 1;

        /* Check for optimality of current solution. */
        descent = -TOLERANCE;
        for (k = 0; k < 4; k++)
        {
            /* Need to find the most negative Lagrange multiplier. */
            if (lamda[k] < descent)
            {
                dIndex = k;
                descent = lamda[k];

                /* Not optimal. */
                optimal = 0;
            }
        }

        if (optimal == 1)
        {
        
			/* Solution is optimal.  Compute dInv if requested. */
			if (dInvReq == 1)
			{
                /* Compute search direction "dInv". */
                cf1000 = gammaVals[1][0] - gammaVals[0][0];
                cf2000 = gammaVals[2][0] - gammaVals[0][0];
                cf2010 = gammaVals[2][0] - gammaVals[1][0];
                cf3000 = gammaVals[3][0] - gammaVals[0][0];
                cf3010 = gammaVals[3][0] - gammaVals[1][0];
                cf3020 = gammaVals[3][0] - gammaVals[2][0];

                dInv[0][0] = (cfArray[1][2] - cfArray[3][2] -
                              cfArray[1][3]) * overDet;
                dInv[0][1] = (gammaVals[2][2] * cf3010 -
                              gammaVals[1][2] * cf3020 -
                              gammaVals[3][2] * cf2010) * overDet;
                dInv[0][2] = (gammaVals[3][1] * cf2010 -
                              gammaVals[2][1] * cf3010 +
                              gammaVals[1][1] * cf3020) * overDet;
				dInv[0][3] = -lamda[0];
                dInv[1][0] = (cfArray[0][3] - cfArray[0][2] +
                              cfArray[3][2]) * overDet;
                dInv[1][1] = (gammaVals[3][2] * cf2000 -
                              gammaVals[2][2] * cf3000 +
                              gammaVals[0][2] * cf3020) * overDet;
                dInv[1][2] = (gammaVals[2][1] * cf3000 -
                              gammaVals[0][1] * cf3020 -
                              gammaVals[3][1] * cf2000) * overDet;
				dInv[1][3] = -lamda[1];
                dInv[2][0] = (cfArray[1][3] - cfArray[0][3] +
                              cfArray[0][1]) * overDet;
                dInv[2][1] = (gammaVals[1][2] * cf3000 -
                              gammaVals[0][2] * cf3010 -
                              gammaVals[3][2] * cf1000) * overDet;
                dInv[2][2] = (gammaVals[0][1] * cf3010 -
                              gammaVals[1][1] * cf3000 +
                              gammaVals[3][1] * cf1000) * overDet;
				dInv[2][3] = -lamda[2];
                dInv[3][0] = (cfArray[0][2] - cfArray[0][1] -
                              cfArray[1][2]) * overDet;
                dInv[3][1] = (gammaVals[2][2] * cf1000 -
                              gammaVals[1][2] * cf2000 +
                              gammaVals[0][2] * cf2010) * overDet;
                dInv[3][2] = (gammaVals[1][1] * cf2000 -
                              gammaVals[0][1] * cf2010 -
                              gammaVals[2][1] * cf1000) * overDet;
				dInv[3][3] = -lamda[3];
	        }

            return 1;
        }

        /* Since there is at least one component of the 
           Lagrange multiplier that is negative, the solution
           is not optimal.  Therefore, we drop the corresponding
           face and try to pick up another.  First, lets find the
           search direction "deltaX" due to the remaining 3
           constraints. */
        switch (dIndex)
        {
            case 0:
                /* Compute search direction "deltaX". */
                deltaX[0] = (cfArray[1][2] - cfArray[3][2] -
                             cfArray[1][3]) * overDet;
                cf1 = gammaVals[3][0] - gammaVals[2][0];
                cf2 = gammaVals[3][0] - gammaVals[1][0];
                cf3 = gammaVals[2][0] - gammaVals[1][0];
                deltaX[1] = (gammaVals[2][2] * cf2 -
                             gammaVals[1][2] * cf1 -
                             gammaVals[3][2] * cf3) * overDet;
                deltaX[2] = (gammaVals[3][1] * cf3 -
                             gammaVals[2][1] * cf2 +
                             gammaVals[1][1] * cf1) * overDet;

                break;

            case 1:
                /* Compute search direction "deltaX". */
                deltaX[0] = (cfArray[0][3] - cfArray[0][2] +
                             cfArray[3][2]) * overDet;
                cf1 = gammaVals[3][0] - gammaVals[2][0];
                cf2 = gammaVals[3][0] - gammaVals[0][0];
                cf3 = gammaVals[2][0] - gammaVals[0][0];
                deltaX[1] = (gammaVals[3][2] * cf3 -
                             gammaVals[2][2] * cf2 +
                             gammaVals[0][2] * cf1) * overDet;
                deltaX[2] = (gammaVals[2][1] * cf2 -
                             gammaVals[0][1] * cf1 -
                             gammaVals[3][1] * cf3) * overDet;

                break;

            case 2:
                /* Compute search direction "deltaX". */
                deltaX[0] = (cfArray[1][3] - cfArray[0][3] +
                             cfArray[0][1]) * overDet;
                cf1 = gammaVals[3][0] - gammaVals[1][0];
                cf2 = gammaVals[3][0] - gammaVals[0][0];
                cf3 = gammaVals[1][0] - gammaVals[0][0];
                deltaX[1] = (gammaVals[1][2] * cf2 -
                             gammaVals[0][2] * cf1 -
                             gammaVals[3][2] * cf3) * overDet;
                deltaX[2] = (gammaVals[0][1] * cf1 -
                             gammaVals[1][1] * cf2 +
                             gammaVals[3][1] * cf3) * overDet;

                break;

            case 3:
                /* Compute search direction "deltaX". */
                deltaX[0] = (cfArray[0][2] - cfArray[0][1] -
                             cfArray[1][2]) * overDet;
                cf1 = gammaVals[2][0] - gammaVals[1][0];
                cf2 = gammaVals[2][0] - gammaVals[0][0];
                cf3 = gammaVals[1][0] - gammaVals[0][0];
                deltaX[1] = (gammaVals[2][2] * cf3 -
                             gammaVals[1][2] * cf2 +
                             gammaVals[0][2] * cf1) * overDet;
                deltaX[2] = (gammaVals[1][1] * cf2 -
                             gammaVals[0][1] * cf1 -
                             gammaVals[2][1] * cf3) * overDet;
        }

        /* Store (x    - P ).
                   init   A
        */
        xMinusPA[0] = result[0] - p1[0];
        xMinusPA[1] = result[1] - p1[1];
        xMinusPA[2] = result[2] - p1[2];

        /* Store (x    - P ).
                   init   B
        */
        xMinusPB[0] = result[0] - p2[0];
        xMinusPB[1] = result[1] - p2[1];
        xMinusPB[2] = result[2] - p2[2];


        /* Trap condition whereby invalid face is dropped, resulting
           in a set of vertices from one object only. */
        if ((*nA == 1) && (jIndex[dIndex][1] == 1))
        {
            
	        /*rotate the deltaX and xMinusP to the original coordinate .*/
            get_origin_vector(deltaX, origin_deltaX ,r1);
            get_origin_vector(xMinusPA, origin_xMinusPA, r1); 

		    /* Compute index into adjacent list of active face from
               object 1. */
            adjPtr = adjSurf1 + jIndex[dIndex][0] * adjCol1;

            /* We start checking adjacent faces of active face of
               object 1.  Set "alpha" to a ridiculously small value. */
            reference = alpha = -1e10;

            /* Make sure we know where we got our "alpha"
               value. */
            whichIndex = 1;

            /* This version does not take into account that some
               adjacent faces may overlap.  It simply make sure
               that none of the faces tested belongs to the active
               set. */
            do
            {
                /* Skip face that was dropped. */
                if (*adjPtr == jIndex[dIndex][0])
                {
                    /* Go on to next face. */
                    *adjPtr++;

                    continue;
                }

                /*                                    -1
                   Compute "bottom" = Gamma .dx  . */
                Gamma_Ptr = Gamma1 + 3*(*adjPtr);
			    returnVal[0] = * Gamma_Ptr ++;
			    returnVal[1] = * Gamma_Ptr ++;
			    returnVal[2] = * Gamma_Ptr ;
                
                top = returnVal[0] * origin_deltaX[0] +
                      returnVal[1] * origin_deltaX[1] +
                      returnVal[2] * origin_deltaX[2] + descent;

                /*           T
                   Compute (a .x     - b )
                             i  init    i
                                 T
                         = (gamma .(x     - P ) - sigma)
                                     init    A
                */
                bottom = returnVal[0] * origin_xMinusPA[0] +
                         returnVal[1] * origin_xMinusPA[1] +
                         returnVal[2] * origin_xMinusPA[2] - result[3];
                value = top / bottom;

                if (value > reference)
                {
                    alpha = value;
                    surface = *adjPtr;
                    reference = alpha ;

                    /* Compute new index into "adjSurf". */
                    adjPtr = adjSurf1 + surface * adjCol1;
                }
                else
                {
                    /* Next adjacent surface. */
                    adjPtr++;
                }
            } while (*adjPtr != -1);

            /* "Alpha" value must be inverted. */
            alpha = 1.0 / alpha;
        }
        else if ((*nB == 1) && (jIndex[dIndex][1] == 2))
        {
            
			/*rotate the deltaX and xMinusP to the original coordinate .*/
            get_origin_vector(deltaX, origin_deltaX ,r2);
            get_origin_vector(xMinusPB, origin_xMinusPB, r2); 
		
		    /* Compute index into adjacent list of active face from
               object 2. */
            adjPtr = adjSurf2 + jIndex[dIndex][0] * adjCol2;
        
            /* Set "alpha" to a ridiculously small value. */
            reference = alpha = -1e10;

            /* Make sure we know where we got our "alpha"
               value. */
            whichIndex = 2;

            /* This version does not take into account that some
               adjacent surfaces may overlap.  It just make sure
               that none of the faces tested belongs to the active
               set. */
            do
            {
                /* Skip face that was dropped. */
                if (*adjPtr == jIndex[dIndex][0])
                {
                    /* Go on to next face. */
                    *adjPtr++;

                    continue;
                }

                /*                                    -1
               Compute "bottom" = Gamma .dx  . */
               Gamma_Ptr = Gamma2 + 3*(*adjPtr);
	   		   returnVal[0] = * Gamma_Ptr ++;
			   returnVal[1] = * Gamma_Ptr ++;
			   returnVal[2] = * Gamma_Ptr ;
	
               top =  returnVal[0] * origin_deltaX[0] +
                      returnVal[1] * origin_deltaX[1] +
                      returnVal[2] * origin_deltaX[2] + descent;

                /*           T
                   Compute (a .x     - b ).
                             i  init    i
                         = (gamma .(x     - P ) - sigma)
                                     init    B
                */
                bottom = returnVal[0] * origin_xMinusPB[0] +
                         returnVal[1] * origin_xMinusPB[1] +
                         returnVal[2] * origin_xMinusPB[2] - result[3];
                value = top / bottom;

                if (value > reference)
                {
                    alpha = value;
                    surface = *adjPtr;
                    reference = alpha ;

                    /* Compute new index into "adjSurf". */
                    adjPtr = adjSurf2 + surface * adjCol2;
                }
                else
                {
                    /* Next adjacent surface. */
                    adjPtr++;
                }
            } while (*adjPtr != -1);
        }
        else
        {
            /* Find out which two faces to check. */
            obj1ActiveFaceCount = obj2ActiveFaceCount = 0;
            obj1AdjFaceCount = obj2AdjFaceCount = INT_MAX;
            for (i = 0; i < 4; i++)
            {
                if (i == dIndex)
                    continue;
                /*To find the active face for each object, and
				these two facea have fewer nieghboring faces.*/
                if (jIndex[i][1] == 1)
                {
                    if (adjFaceCount1[jIndex[i][0]] < obj1AdjFaceCount)
                    {
                        obj1AdjFaceCount = adjFaceCount1[jIndex[i][0]];
                        obj1Face = jIndex[i][0];
                        ++obj1ActiveFaceCount;
                    }
                }
                else
                {
                    if (adjFaceCount2[jIndex[i][0]] < obj2AdjFaceCount)
                    {
                        obj2AdjFaceCount = adjFaceCount2[jIndex[i][0]];
                        obj2Face = jIndex[i][0];
                        ++obj2ActiveFaceCount;
                    }
                }
            }

            if (obj1ActiveFaceCount == 1)
            {
                /* Compute index into adjacent list of active face from
                   object 1. */
                adjPtr = adjSurf1 + obj1Face * adjCol1;

                /* We start checking adjacent faces of active face of
                   object 1. */
                if (obj1AdjFaceCount > THRESHOLD)
                {
                    binarySearchFinal(getGamma1, Gamma1, r1, deltaX,
                                      xMinusPA, result[3], descent, adjPtr,
                                      obj1AdjFaceCount, &surface, &alpha);
                    reference = alpha ;
                    whichIndex = 1;
                }
                else
                {
                    /* Set "alpha" to a ridiculously large value. */
                    reference = alpha = 1e10;

                    /* This version does not take into account that some
                       adjacent faces may overlap.  It simply make sure
                       that none of the faces tested belongs to the active
                       set. */
                    
					/*rotate the deltaX and xMinusP to the original coordinate .*/
                    get_origin_vector(deltaX, origin_deltaX ,r1);
                    get_origin_vector(xMinusPA, origin_xMinusPA, r1); 
					
					do
                    {
                        if (inIndex(*adjPtr, iIndex1, *nA) != 1)
                        {
                            /*                                    -1
                               Compute "bottom" = Gamma .dx  . */
                            Gamma_Ptr = Gamma1 + 3*(*adjPtr);
			                returnVal[0] = * Gamma_Ptr ++;
			                returnVal[1] = * Gamma_Ptr ++;
			                returnVal[2] = * Gamma_Ptr ;

                            bottom = returnVal[0] * origin_deltaX[0] +
                                     returnVal[1] * origin_deltaX[1] +
                                     returnVal[2] * origin_deltaX[2] + descent;


                            if (bottom < 0.0)
                            {
                                /*           T
                                   Compute (a .x     - b )
                                             i  init    i
                                                 T
                                         = (gamma .(x     - P ) - sigma)
                                                     init    A
                                */
                                top = returnVal[0] * xMinusPA[0] +
                                      returnVal[1] * xMinusPA[1] +
                                      returnVal[2] * xMinusPA[2] - result[3];
                                value = top / bottom;

                                if (value < reference)
                                {
                                    alpha = value;
                                    surface = *adjPtr;
                                    reference = alpha ;

                                    /* Make sure we know where we got our "alpha"
                                       value. */
                                    whichIndex = 1;
                                }
                            }
                        }

                        /* Next adjacent surface. */
                        adjPtr++;

                    } while (*adjPtr != -1);
                }

                /* Compute index into adjacent list of active face from
                   object 2. */
                adjPtr = adjSurf2 + obj2Face * adjCol2;
        
                /* This version does not take into account that some
                   adjacent surfaces may overlap.  It just make sure
                   that none of the faces tested belongs to the active
                   set. */
	            /*rotate the deltaX and xMinusP to the original coordinate .*/
                get_origin_vector(deltaX, origin_deltaX ,r2);
                get_origin_vector(xMinusPB, origin_xMinusPB, r2); 

                do
                {
                    if (inIndex(*adjPtr, iIndex2, *nB) != 1)
                    {
	                /*                                    -1
                           Compute "bottom" = Gamma .dx  . */
                        Gamma_Ptr = Gamma2 + 3*(*adjPtr);
			            returnVal[0] = * Gamma_Ptr ++;
			            returnVal[1] = * Gamma_Ptr ++;
			            returnVal[2] = * Gamma_Ptr ;

                        bottom = returnVal[0] * origin_deltaX[0] +
                                 returnVal[1] * origin_deltaX[1] +
                                 returnVal[2] * origin_deltaX[2] + descent;

                        if (bottom < 0.0)
                        {
                            /*           T
                               Compute (a .x     - b ).
                                         i  init    i
                                     = (gamma .(x     - P ) - sigma)
                                                 init    B
                            */
                            top = returnVal[0] * origin_xMinusPB[0] +
                                  returnVal[1] * origin_xMinusPB[1] +
                                  returnVal[2] * origin_xMinusPB[2] - result[3];
                            value = top / bottom;

                            if (value < reference)
                            {
                                alpha = value;
                                surface = *adjPtr;
                                reference = alpha ;

                                /* Make sure we know where we got our "alpha"
                                   value. */
                                whichIndex = 2;
                            }
                        }
                    }

                    /* Next adjacent surface. */
                    adjPtr++;

                } while (*adjPtr != -1);
            }
            else
            {
                /* Compute index into adjacent list of active face from
                   object 2. */
                adjPtr = adjSurf2 + obj2Face * adjCol2;

                if (obj2AdjFaceCount > THRESHOLD)
                {
                    binarySearchFinal(getGamma2, Gamma2, r2, deltaX,
                                      xMinusPB, result[3], descent, adjPtr,
                                      obj2AdjFaceCount, &surface, &alpha);
                    reference = alpha ;
                    whichIndex = 2;
                }
                else
                {
                    /* Set "alpha" to a ridiculously large value. */
                    reference = alpha = 1e10;

                    /* This version does not take into account that some
                       adjacent surfaces may overlap.  It just make sure
                       that none of the faces tested belongs to the active
                       set. */
                    /*rotate the deltaX and xMinusP to the original coordinate .*/
                    get_origin_vector(deltaX, origin_deltaX ,r2);
                    get_origin_vector(xMinusPB, origin_xMinusPB, r2); 

                    do
                    {
                        if (inIndex(*adjPtr, iIndex2, *nB) != 1)
                        {
	                    /*                                    -1
                               Compute "bottom" = Gamma .dx  . */
                            Gamma_Ptr = Gamma2 + 3*(*adjPtr);
			                returnVal[0] = * Gamma_Ptr ++;
			                returnVal[1] = * Gamma_Ptr ++;
			                returnVal[2] = * Gamma_Ptr ;

                            bottom = returnVal[0] * origin_deltaX[0] +
                                     returnVal[1] * origin_deltaX[1] +
                                     returnVal[2] * origin_deltaX[2] + descent;

                            if (bottom < 0.0)
                            {
                                /*           T
                                   Compute (a .x     - b ).
                                             i  init    i
                                         = (gamma .(x     - P ) - sigma)
                                                     init    B
                                */
                                top = returnVal[0] * origin_xMinusPB[0] +
                                      returnVal[1] * origin_xMinusPB[1] +
                                      returnVal[2] * origin_xMinusPB[2] - result[3];
                                value = top / bottom;

                                if (value < reference)
                                {
                                    alpha = value;
                                    surface = *adjPtr;
                                    reference = alpha ;

                                    /* Make sure we know where we got our "alpha"
                                       value. */
                                    whichIndex = 2;
                                }
                            }
                        }

                        /* Next adjacent surface. */
                        adjPtr++;

                    } while (*adjPtr != -1);
                }

                /* Compute index into adjacent list of active face from
                   object 1. */
                adjPtr = adjSurf1 + obj1Face * adjCol1;

                /* This version does not take into account that some
                   adjacent faces may overlap.  It simply make sure
                   that none of the faces tested belongs to the active
                   set. */
	            /*rotate the deltaX and xMinusP to the original coordinate .*/
                get_origin_vector(deltaX, origin_deltaX ,r1);
                get_origin_vector(xMinusPA, origin_xMinusPA, r1); 

                do
                {
                    if (inIndex(*adjPtr, iIndex1, *nA) != 1)
                    {

                        /*                                    -1
                           Compute "bottom" = Gamma .dx  . */
                        Gamma_Ptr = Gamma1 + 3*(*adjPtr);
			            returnVal[0] = * Gamma_Ptr ++;
			            returnVal[1] = * Gamma_Ptr ++;
			            returnVal[2] = * Gamma_Ptr ;

                        bottom = returnVal[0] * origin_deltaX[0] +
                                 returnVal[1] * origin_deltaX[1] +
                                 returnVal[2] * origin_deltaX[2] + descent;

                        if (bottom < 0.0)
                        {
                            /*           T
                               Compute (a .x     - b )
                                         i  init    i
                                             T
                                     = (gamma .(x     - P ) - sigma)
                                                 init    A
                            */
                            top = returnVal[0] * origin_xMinusPA[0] +
                                  returnVal[1] * origin_xMinusPA[1] +
                                  returnVal[2] * origin_xMinusPA[2] - result[3];
                            value = top / bottom;

                            if (value < reference)
                            {
                                alpha = value;
                                surface = *adjPtr;
                                reference = alpha ;

                                /* Make sure we know where we got our "alpha"
                                   value. */
                                whichIndex = 1;
                            }
                        }
                    }

                    /* Next adjacent surface. */
                    adjPtr++;

                } while (*adjPtr != -1);

            }
        }

        /* Update "iIndex". */
        if (jIndex[dIndex][1] == whichIndex)
        {
            /* The face dropped and the new one picked up are from
               the same object.  We need to compact "iIndex" and insert
               the new value in front. */
            if (whichIndex == 1)
            {
                /* Search for occurrence of surface to be replaced
                   in "iIndex1". */
                switch (*nA)
                {
                    case 3:
                        if (iIndex1[2] == jIndex[dIndex][0])
                        {
                            iIndex1[2] = iIndex1[1];
                            iIndex1[1] = iIndex1[0];
                            iIndex1[0] = surface;

                            break;
                        }
                    case 2:
                        if (iIndex1[1] == jIndex[dIndex][0])
                        {
                            iIndex1[1] = iIndex1[0];
                            iIndex1[0] = surface;

                            break;
                        }
                    default:
                        iIndex1[0] = surface;
                }
            }
            else
            {
                /* Search for occurrence of surface to be replaced
                   in "iIndex2". */
                switch (*nB)
                {
                    case 3:
                        if (iIndex2[2] == jIndex[dIndex][0])
                        {
                            iIndex2[2] = iIndex2[1];
                            iIndex2[1] = iIndex2[0];
                            iIndex2[0] = surface;

                            break;
                        }
                    case 2:
                        if (iIndex2[1] == jIndex[dIndex][0])
                        {
                            iIndex2[1] = iIndex2[0];
                            iIndex2[0] = surface;

                            break;
                        }
                    default:
                        iIndex2[0] = surface;
                }
            }
        }
        else
        {
            /* The face dropped and the new one picked up are from
               different objects.  More work need to be done here. */
            if (whichIndex == 1)
            {
                /* New surface is from object 1. */
                switch (*nA)
                {
                    case 2: iIndex1[2] = iIndex1[1];
                    case 1: iIndex1[1] = iIndex1[0];
                            iIndex1[0] = surface;
                }
                (*nA)++;

                /* Remove old surface from "iIndex2". */
                for (i = 0; i < (*nB - 1); i++)
                {
                    if (iIndex2[i] == jIndex[dIndex][0])
                    {
                        switch (i)
                        {
                            case 0: iIndex2[0] = iIndex2[1];
                            case 1: iIndex2[1] = iIndex2[2];
                        }

                        break;
                    }
                }
                (*nB)--;
            }
            else
            {
                /* New surface is from object 2. */
                switch (*nB)
                {
                    case 2: iIndex2[2] = iIndex2[1];
                    case 1: iIndex2[1] = iIndex2[0];
                            iIndex2[0] = surface;
                }
                (*nB)++;

                /* Remove old surface from "iIndex1". */
                for (i = 0; i < (*nA - 1); i++)
                {
                    if (iIndex1[i] == jIndex[dIndex][0])
                    {
                        switch (i)
                        {
                            case 0: iIndex1[0] = iIndex1[1];
                            case 1: iIndex1[1] = iIndex1[2];
                        }

                        break;
                    }
                }
                (*nA)--;
            }
        }

        /* Update "jIndex". */
        jIndex[dIndex][0] = surface;
        jIndex[dIndex][1] = whichIndex;

        /* Compute result of run. */
        result[0] -= alpha * deltaX[0];
        result[1] -= alpha * deltaX[1];
        result[2] -= alpha * deltaX[2];
        result[3] += alpha * descent;

        /* Update stored gamma values. */
        if (whichIndex == 1)
            getGamma1(Gamma1, r1, surface, gammaVals[dIndex]);
        else
            getGamma2(Gamma2, r2, surface, gammaVals[dIndex]);

        /* Recompute all Lagrange multipliers. */
        switch (dIndex)
        {
            case 0:
                /* Need to recompute "det012", "det013" and "det032". */
                cfArray[0][1] = gammaVals[0][1] * gammaVals[1][2] -
                                gammaVals[0][2] * gammaVals[1][1];
                cfArray[0][2] = gammaVals[0][1] * gammaVals[2][2] -
                                gammaVals[0][2] * gammaVals[2][1];
                cfArray[0][3] = gammaVals[0][1] * gammaVals[3][2] -
                                gammaVals[0][2] * gammaVals[3][1];

                det012 = gammaVals[0][0] * cfArray[1][2] -
                         gammaVals[1][0] * cfArray[0][2] +
                         gammaVals[2][0] * cfArray[0][1];
                det013 = gammaVals[0][0] * cfArray[1][3] -
                         gammaVals[1][0] * cfArray[0][3] +
                         gammaVals[3][0] * cfArray[0][1];
                det032 = gammaVals[0][0] * cfArray[3][2] -
                         gammaVals[3][0] * cfArray[0][2] +
                         gammaVals[2][0] * cfArray[0][3];

                break;

            case 1:
                /* Need to recompute "det012", "det013" and "det312". */
                cfArray[0][1] = gammaVals[0][1] * gammaVals[1][2] -
                                gammaVals[0][2] * gammaVals[1][1];
                cfArray[1][2] = gammaVals[1][1] * gammaVals[2][2] -
                                gammaVals[1][2] * gammaVals[2][1];
                cfArray[1][3] = gammaVals[1][1] * gammaVals[3][2] -
                                gammaVals[1][2] * gammaVals[3][1];

                det012 = gammaVals[0][0] * cfArray[1][2] -
                         gammaVals[1][0] * cfArray[0][2] +
                         gammaVals[2][0] * cfArray[0][1];
                det013 = gammaVals[0][0] * cfArray[1][3] -
                         gammaVals[1][0] * cfArray[0][3] +
                         gammaVals[3][0] * cfArray[0][1];
                det312 = gammaVals[3][0] * cfArray[1][2] -
                         gammaVals[1][0] * cfArray[3][2] -
                         gammaVals[2][0] * cfArray[1][3];

                break;

            case 2:
                /* Need to recompute "det012", "det032" and "det312". */
                cfArray[0][2] = gammaVals[0][1] * gammaVals[2][2] -
                                gammaVals[0][2] * gammaVals[2][1];
                cfArray[1][2] = gammaVals[1][1] * gammaVals[2][2] -
                                gammaVals[1][2] * gammaVals[2][1];
                cfArray[3][2] = gammaVals[3][1] * gammaVals[2][2] -
                                gammaVals[3][2] * gammaVals[2][1];

                det012 = gammaVals[0][0] * cfArray[1][2] -
                         gammaVals[1][0] * cfArray[0][2] +
                         gammaVals[2][0] * cfArray[0][1];
                det032 = gammaVals[0][0] * cfArray[3][2] -
                         gammaVals[3][0] * cfArray[0][2] +
                         gammaVals[2][0] * cfArray[0][3];
                det312 = gammaVals[3][0] * cfArray[1][2] -
                         gammaVals[1][0] * cfArray[3][2] -
                         gammaVals[2][0] * cfArray[1][3];

                break;

            case 3:
                /* Need to recompute "det013", "det032" and "det312". */
                cfArray[0][3] = gammaVals[0][1] * gammaVals[3][2] -
                                gammaVals[0][2] * gammaVals[3][1];
                cfArray[1][3] = gammaVals[1][1] * gammaVals[3][2] -
                                gammaVals[1][2] * gammaVals[3][1];
                cfArray[3][2] = gammaVals[3][1] * gammaVals[2][2] -
                                gammaVals[3][2] * gammaVals[2][1];

                det013 = gammaVals[0][0] * cfArray[1][3] -
                         gammaVals[1][0] * cfArray[0][3] +
                         gammaVals[3][0] * cfArray[0][1];
                det032 = gammaVals[0][0] * cfArray[3][2] -
                         gammaVals[3][0] * cfArray[0][2] +
                         gammaVals[2][0] * cfArray[0][3];
                det312 = gammaVals[3][0] * cfArray[1][2] -
                         gammaVals[1][0] * cfArray[3][2] -
                         gammaVals[2][0] * cfArray[1][3];
        }

           /*calculate test condition.*/
	    Max_det = fabs(det012);
	    if ( fabs(det013) > Max_det) Max_det = fabs(det013);
	    if ( fabs(det032) > Max_det) Max_det = fabs(det032);	
	    if ( fabs(det312) > Max_det) Max_det = fabs(det312);
        testCondition = fabs(bottom)/Max_det;

        /* Check for divide-by-zero. */
        bottom = det012 - det013 - det032 - det312;
        if (testCondition < TOLERANCE)
            return 3;

        overDet = 1.0 / bottom;

        lamda[0] = -(det312 * overDet);
        lamda[1] = -(det032 * overDet);
        lamda[2] = -(det013 * overDet);
        lamda[3] = 1.0 - (lamda[0] + lamda[1] + lamda[2]);

    }
}



/* ================================================================ */
/* Function to check if the listed surface is stored in the given
   "iIndex" array.


   Input Parameters
   ================
   indexValue -- Surface number to check for.
   *iIndexPtr -- Pointer to adjacent surface list to be checked.
   size       -- Number of elements in "iIndexPtr" array.


   Return Value
   ============
   The function returns a 1 if the surface number given by "indexValue"
   is found in the index "iIndexPtr", otherwise a value of 0 is returned.
 */
int inIndex(int indexValue, int iIndexPtr[], int size)
{
    int *i;

    for (i = iIndexPtr; i < (iIndexPtr + size); i++)
    {
        if (*i == indexValue)
            /* Found the surface. */
            return 1;
    }

    /* No surface found. */
    return 0;
}



/* ================================================================ */
/*                                                   T
   Function to compute the solution of the equation A x = b, where
    T                                      T
   A  is a 3 x 3 matrix and b = (-1 -1 -1) .


   Input Parameters
   ================
   A[3][3]     -- The 3 x 3 matrix A.  Note, although the
                  matrix A is passed in, this routine obtains the
                               T
                  solution of A x = b, not Ax = b.
   i0, i1, i2  -- Inidices indicating which three columns of matrix
                  A is to be used for solving the linear equation.


   Output Parameter
   ================
   x[3]        -- 3 x 1 vector to hold the result.


   Return Parameter
   ================
   A 0 indicates a solution is obtained.  Otherwise a 1 is returned
   to indicate that the matrix A given is singular.
*/
__inline int linearSolver3x3c(double A[][3], int i0, int i1,
                              int i2, double x[])
{
    /* Variables to store computed determinants. */
    double det012, det013, det032, det312;

    /* Array to hold computed cofactors so as to eliminate
       redundant computations. */
    double cfArray[4][4];

    /* Miscellaneous variables. */
    double overDet;


    cfArray[0][1] = A[i1][0] * A[i2][1] - A[i2][0] * A[i1][1];
    cfArray[0][2] = A[i1][0] * A[i2][2] - A[i2][0] * A[i1][2];

    /*                               T
       Remember, since b = (-1 -1 -1) , any cofactors that involves
       "b" needs no real multiplications. */
    cfArray[0][3] = A[i2][0] - A[i1][0];
    cfArray[1][2] = A[i1][1] * A[i2][2] - A[i2][1] * A[i1][2];
    cfArray[1][3] = A[i2][1] - A[i1][1];
    cfArray[3][2] = A[i1][2] - A[i2][2];

    det012 = A[i0][0] * cfArray[1][2] - A[i0][1] * cfArray[0][2] +
             A[i0][2] * cfArray[0][1];

    if ((det012 > -TOLERANCE) && (det012 < TOLERANCE))
        /* Singular A. */
        return 1;

    overDet = 1.0 / det012;

    det013 = A[i0][0] * cfArray[1][3] - A[i0][1] * cfArray[0][3] -
             cfArray[0][1];
    det032 = A[i0][0] * cfArray[3][2] + cfArray[0][2] +
             A[i0][2] * cfArray[0][3];
    det312 = -(cfArray[1][2] + A[i0][1] * cfArray[3][2] +
               A[i0][2] * cfArray[1][3]);

    x[0] = det312 * overDet;
    x[1] = det032 * overDet;
    x[2] = det013 * overDet;

    return 0;
}

/* ================================================================ */
/* Function to get the original coordinate of a vector, the current rotation matrix                                                   
is r .*/
void get_origin_vector(double x[], double origin_x[], double r[][3])
{
	origin_x[0] = r[0][0]*x[0]+r[1][0]*x[1]+r[2][0]*x[2];
	origin_x[1] = r[0][1]*x[0]+r[1][1]*x[1]+r[2][1]*x[2];
	origin_x[2] = r[0][2]*x[0]+r[1][2]*x[1]+r[2][2]*x[2];
}