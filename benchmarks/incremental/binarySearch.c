#include "Gamma.h"
#include "binarySearch.h"
#include "Solver4c.h"


/* Definitions used to indicate gradient. */
#define POSITIVE 0
#define NEGATIVE 1


/* Local prototype declaration. */
__inline double getMuInv(void (*getGamma)(), double *Gamma,
                         double origin_deltaX[],
                         double origin_xMinusP[], double sigma, int face);
__inline double getMuInvFinal(void (*getGamma)(),double *Gamma,
							  double origin_deltaX[], double origin_xMinusP[],
                              double sigma, double descent, int face);


/* ================================================================ */
/* This function implements a pseudo-binary search to obtain the
   maximum feasible step size for further descent. */
int binarySearch(void (*getGamma)(), double *Gamma,
                 double rotate[][3], double deltaX[],
                 double xMinusP[], double sigma, int *adjSurf,
                 int adjFaceCount, int *activeFace,
                 double *stepSize)
{
    /* Control variables. */
    int gradLowBound, gradTrialPoint;
    int lowBound, upBound, trialPoint;
    double lowBoundVal, upBoundVal, trialPointVal;
	double origin_xMinusP[3], origin_deltaX[3];
	

    /* Miscellaneous variables. */
    double t1, t2;

	/*rotate the deltaX and xMinusP to the original coordinate .*/
    get_origin_vector(deltaX, origin_deltaX, rotate);
    get_origin_vector(xMinusP, origin_xMinusP, rotate); 

    
	/* ====================== */
    /* Setup initial zone that covers the complete face. */
    lowBound = 0;
    lowBoundVal = getMuInv(getGamma, Gamma,
		                   origin_deltaX, origin_xMinusP,
                           sigma, adjSurf[lowBound]);
    
	upBound = adjFaceCount - 1;
    upBoundVal = getMuInv(getGamma, Gamma,
		                  origin_deltaX, origin_xMinusP,
                          sigma, adjSurf[upBound]);

    /* First make sure that "lowBound" is not a maxima. */
    t1 = upBoundVal;
    t2 = getMuInv(getGamma, Gamma,
		          origin_deltaX, origin_xMinusP,
                  sigma, adjSurf[lowBound + 1]);
    if ((t1 <= lowBoundVal) && (t2 <= lowBoundVal))
    {
        /* Solution is optimal. */
        *activeFace = adjSurf[lowBound];
        *stepSize = 1.0 / lowBoundVal;

        return 0;
    }
    
    /* Obtain gradient of "lowBound". */
    if (lowBoundVal < t2)
        gradLowBound = POSITIVE;
    else
        gradLowBound = NEGATIVE;

    /* Check that "upBound" is not a maxima too. */
    t1 = getMuInv(getGamma, Gamma,origin_deltaX, origin_xMinusP,
                  sigma, adjSurf[upBound - 1]);
    t2 = lowBoundVal;
    if ((t1 <= upBoundVal) && (t2 <= upBoundVal))
    {
        /* Solution is optimal. */
        *activeFace = adjSurf[upBound];
        *stepSize = 1.0 / upBoundVal;

        return 0;
    }

    /* Unwrapped recursive function to do "pseudo binary searches". */
    while (1)
    {
        /* Start with left test zone, i.e. comprising "lowerBound" to
           "trialPoint". */
        trialPoint = (lowBound + upBound) >> 1;
        trialPointVal = getMuInv(getGamma, Gamma,origin_deltaX,
                                 origin_xMinusP, sigma, adjSurf[trialPoint]);

        /* Check that "trialPoint" is not a maxima. */
        t1 = getMuInv(getGamma, Gamma,origin_deltaX, origin_xMinusP,
                      sigma, adjSurf[trialPoint - 1]);
        t2 = getMuInv(getGamma, Gamma,origin_deltaX, origin_xMinusP,
                      sigma, adjSurf[trialPoint + 1]);
        if ((t1 <= trialPointVal) && (t2 <= trialPointVal))
        {
            /* Solution is optimal. */
            *activeFace = adjSurf[trialPoint];
            *stepSize = 1.0 / trialPointVal;
    
            return 0;
        }
        
        /* Obtain gradient of "trialPoint". */
        if (trialPointVal < t2)
            gradTrialPoint = POSITIVE;
        else
            gradTrialPoint = NEGATIVE;

        /* Binary search logic. */
        if (gradLowBound == POSITIVE)
        {
            if (gradTrialPoint == NEGATIVE)
                upBound = trialPoint;
            else
            {
                /* "gradTrialPoint" is positive, so check values of "lowBound"
                   and "trialPoint" to decide which to zone to pick up. */
                if (lowBoundVal < trialPointVal)
                {
                    lowBound = trialPoint;
                    lowBoundVal = trialPointVal;
                    gradLowBound = gradTrialPoint;
                }
                else
                    upBound = trialPoint;
            }
        }
        else
        {
            if (gradTrialPoint == POSITIVE)
            {
                lowBound = trialPoint;
                lowBoundVal = trialPointVal;
                gradLowBound = gradTrialPoint;
            }
            else
            {
                /* "gradTrialPoint" is negative, so check values of "lowBound"
                   and "trialPoint" to decide which to zone to pick up. */
                if (lowBoundVal < trialPointVal)
                    upBound = trialPoint;
                else
                {
                    lowBound = trialPoint;
                    lowBoundVal = trialPointVal;
                    gradLowBound = gradTrialPoint;
                }
            }
        }
    }
}


/* ================================================================ */
/* This function computes the inverse of the "mu" parameter. */
__inline double getMuInv(void (*getGamma)(), 
                         double *Gamma, double origin_deltaX[],
                         double origin_xMinusP[], double sigma, int face)
{
    /* Temporary variables. */
    double bottom, top;

    /* Miscellaneous variables. */
    double returnVal[3];
	double *Gamma_Ptr;


    /*                                       -1
          Compute "bottom" = Gamma .dx  . */
    Gamma_Ptr = Gamma + 3 * face;
	returnVal[0] = * Gamma_Ptr ++;
	returnVal[1] = * Gamma_Ptr ++;
	returnVal[2] = * Gamma_Ptr ;
	
    top = returnVal[0] * origin_deltaX[0] + returnVal[1] *origin_deltaX[1] +
          returnVal[2] * origin_deltaX[2] + 1.0;

    /*                       T
       Compute (sigma - gamma .(x - P)).
    */
    bottom = sigma - (returnVal[0] * origin_xMinusP[0] +
                      returnVal[1] * origin_xMinusP[1] +
                      returnVal[2] * origin_xMinusP[2]);

    /*         1
       Return ----.
               mu
    */
    return (top / bottom);
}



/* ================================================================ */
/* This function implements a pseudo-binary search to obtain the
   maximum feasible step size for further descent in final descent. */
int binarySearchFinal(void (*getGamma)(), double *Gamma,
                      double rotate[][3], double deltaX[], 
                      double xMinusP[], double sigma, double descent,
                      int *adjSurf, int adjFaceCount, int *activeFace,
                      double *stepSize)
{
    /* Control variables. */
    int gradLowBound, gradTrialPoint;
    int lowBound, upBound, trialPoint;
    double lowBoundVal, upBoundVal, trialPointVal;
    double origin_deltaX[3], origin_xMinusP[3];
	

    /* Miscellaneous variables. */
    double t1, t2;

	/*rotate the deltaX and xMinusP to the original coordinate .*/
    get_origin_vector(deltaX, origin_deltaX ,rotate);
    get_origin_vector(xMinusP, origin_xMinusP, rotate); 
	
    /* ====================== */
    /* Setup initial zone that covers the complete face. */
    lowBound = 0;
    lowBoundVal = getMuInvFinal(getGamma, Gamma,
                                origin_deltaX, origin_xMinusP, sigma, descent,
                                adjSurf[lowBound]);
    upBound = adjFaceCount - 1;
    upBoundVal = getMuInvFinal(getGamma, Gamma,origin_deltaX,
                               origin_xMinusP, sigma, descent, adjSurf[upBound]);

    /* First make sure that "lowBound" is not a maxima. */
    t1 = upBoundVal;
    t2 = getMuInvFinal(getGamma, Gamma,origin_deltaX, origin_xMinusP,
                       sigma, descent, adjSurf[lowBound + 1]);
    if ((t1 <= lowBoundVal) && (t2 <= lowBoundVal))
    {
        /* Solution is optimal. */
        *activeFace = adjSurf[lowBound];
        *stepSize = 1.0 / lowBoundVal;

        return 0;
    }
    
    /* Obtain gradient of "lowBound". */
    if (lowBoundVal < t2)
        gradLowBound = POSITIVE;
    else
        gradLowBound = NEGATIVE;

    /* Check that "upBound" is not a maxima too. */
    t1 = getMuInvFinal(getGamma,Gamma,origin_deltaX,
                       origin_xMinusP, sigma, descent, adjSurf[upBound - 1]);
    t2 = lowBoundVal;
    if ((t1 <= upBoundVal) && (t2 <= upBoundVal))
    {
        /* Solution is optimal. */
        *activeFace = adjSurf[upBound];
        *stepSize = 1.0 / upBoundVal;

        return 0;
    }

    /* Unwrapped recursive function to do "pseudo binary searches". */
    while (1)
    {
        /* Start with left test zone, i.e. comprising "lowerBound" to
           "trialPoint". */
        trialPoint = (lowBound + upBound) >> 1;
        trialPointVal = getMuInvFinal(getGamma, Gamma,origin_deltaX,
			                          origin_xMinusP, sigma, descent,
                                      adjSurf[trialPoint]);

        /* Check that "trialPoint" is not a maxima. */
        t1 = getMuInvFinal(getGamma, Gamma,origin_deltaX,
                           origin_xMinusP, sigma, descent,
                           adjSurf[trialPoint - 1]);
        t2 = getMuInvFinal(getGamma, Gamma,origin_deltaX,
                           origin_xMinusP, sigma, descent,
                           adjSurf[trialPoint + 1]);
        if ((t1 <= trialPointVal) && (t2 <= trialPointVal))
        {
            /* Solution is optimal. */
            *activeFace = adjSurf[trialPoint];
            *stepSize = 1.0 / trialPointVal;
    
            return 0;
        }
        
        /* Obtain gradient of "trialPoint". */
        if (trialPointVal < t2)
            gradTrialPoint = POSITIVE;
        else
            gradTrialPoint = NEGATIVE;

        /* Binary search logic. */
        if (gradLowBound == POSITIVE)
        {
            if (gradTrialPoint == NEGATIVE)
                upBound = trialPoint;
            else
            {
                /* "gradTrialPoint" is positive, so check values of "lowBound"
                   and "trialPoint" to decide which to zone to pick up. */
                if (lowBoundVal < trialPointVal)
                {
                    lowBound = trialPoint;
                    lowBoundVal = trialPointVal;
                    gradLowBound = gradTrialPoint;
                }
                else
                    upBound = trialPoint;
            }
        }
        else
        {
            if (gradTrialPoint == POSITIVE)
            {
                lowBound = trialPoint;
                lowBoundVal = trialPointVal;
                gradLowBound = gradTrialPoint;
            }
            else
            {
                /* "gradTrialPoint" is negative, so check values of "lowBound"
                   and "trialPoint" to decide which to zone to pick up. */
                if (lowBoundVal < trialPointVal)
                    upBound = trialPoint;
                else
                {
                    lowBound = trialPoint;
                    lowBoundVal = trialPointVal;
                    gradLowBound = gradTrialPoint;
                }
            }
        }
    }
}


/* ================================================================ */
/* This function computes the inverse of the "mu" parameter. */
__inline double getMuInvFinal(void (*getGamma)(), double *Gamma,
                              double origin_deltaX[], double origin_xMinusP[],
                              double sigma, double descent, int face)
{
    /* Temporary variables. */
    double bottom, top;

    /* Miscellaneous variables. */
    double returnVal[3];
	double *Gamma_Ptr;


    /*                                       -1
          Compute "bottom" = Gamma .dx  . */
    Gamma_Ptr = Gamma + 3 * face;
	returnVal[0] = * Gamma_Ptr ++;
	returnVal[1] = * Gamma_Ptr ++;
	returnVal[2] = * Gamma_Ptr ;
	
    top = returnVal[0] * origin_deltaX[0] + returnVal[1] *origin_deltaX[1] +
          returnVal[2] * origin_deltaX[2] + descent;

   
    /*           T
       Compute (a .x - b )
                     T
             = (gamma .(x - P ) - sigma)
    */
    bottom = returnVal[0] * origin_xMinusP[0] +
             returnVal[1] * origin_xMinusP[1] +
             returnVal[2] * origin_xMinusP[2] - sigma;

    /*         1
       Return ----.
               mu
    */
    return (top / bottom);
}

