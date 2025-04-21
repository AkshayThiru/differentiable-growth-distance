#include <malloc.h>
#include <limits.h>
#include "Gamma.h"



/* ID used to indicate an initialized array. */
#define INIT_ID         0
#define FIRST_RUN_ID    1

/* ID value used by "gGammaIndex" and "gGammaPIndex" to indicate a useful
   value in "gGamma" and "gGammaP".  The use of this ID preempts the
   requirement of reinitializing "gGamma" array to zero. */
int *runID;

/* ========================== */
/* Data structure to hold computed values of gamma together with
   additional data for its proper management. */
struct gammaData
{
    /* Pointer to array holding commputed gamma values. */
    double *gGamma;

    /* Pointer to array holding indicators to indicate if the
       corresponding value stored in "gGamma" is meaningful. */
    int *gGammaIndex;

    /* Size of the above arrays. */
    int size;
};


/* Data structure to hold computed values of gamma together with
   additional data for its proper management. */
struct gammaPData
{
    /*                                          T
       Pointer to array holding commputed (gamma  * p) values.
    */
    double *gGammaP;

    /* Pointer to array holding indicators to indicate if the
       corresponding value stored in "gGammaP" is meaningful. */
    int *gGammaPIndex;

    /* Size of the above arrays. */
    int size;
};



/* ========================== */
/* Global Declarations. */
struct gammaData gGamma1;
struct gammaData gGamma2;
struct gammaPData gGamma1P1, gGamma1P2;
struct gammaPData gGamma2P1, gGamma2P2;


/* ========================== */
/* Local prototypes. */
__inline void getGamma(struct gammaData *gGamma, double *Gamma,
					   double rotate[][3], int index, double result[]);
__inline double getGammaP(struct gammaPData *gGammaP,struct gammaData *gGamma,
                          double *Gamma, double rotate[][3], double p[],
                          int index);

__inline void initgGammaIndex(struct gammaData *gGamma);
__inline void initgGammaPIndex(struct gammaPData *gGammaP);



/* ================================================================ */
/*                                               T
                                                n    T
   This function returns the computed value of ----.R .  To speed
                                                T -
                                               n .x
   up the operation, the computed value is stored in the "gGamma"
   array so that it can be referenced in future.  Further request
   for this computed value results in a simple retrieval operation
   instead of recomputation.


   Input Parameters
   ================
   *gGamma      -- Structure that holds the computed gamma information.
   *normals     -- Matrix containing the set of normals that defines the
                   faces that make up the object.
   *xBars       -- Matrix containing the coordinates of points on the
                   faces of the object corresponding to the normals
                   defined in the matrix "normals".
   rotate[3][3] -- Contains the transpose of the rotation matrix.
   index        -- Index into "normals" and "xBars" matrices for which the
                   corresponding gamma value is required.


   Output Parameter
   ================
   result[3]    -- Requested gamma value.
*/
__inline void getGamma(struct gammaData *gGamma, double *Gamma,
					   double rotate[][3], int index, double result[])
{
    /* Pointers to all required data. */
    double *g_x_Ptr, *g_y_Ptr, *g_z_Ptr;
	double *Gamma_Ptr;
	double origin_x, origin_y, origin_z;

    /* Miscellaneous variables. */
    register int i;

    /* Compute index into array. */
    i = index * 3;

    /* Set "gGamma" pointers. */
    g_x_Ptr = gGamma->gGamma + i;
    g_y_Ptr = g_x_Ptr+1;
	g_z_Ptr = g_y_Ptr+1;
    
    /* Check to see if it is necessary to recompute. */
    if (gGamma->gGammaIndex[index] != (*runID))
    {
        /* Need to compute gamma.  So get the original gamma value first. */
        Gamma_Ptr = Gamma + i;
		origin_x = *Gamma_Ptr++;
	    origin_y = *Gamma_Ptr++;
	    origin_z = *Gamma_Ptr;
        /*           
           Compute  Gamma = rotate[][3].original Gamma.
        */
        /* the rotate[][3] is the transpose of the rotation matrix.*/
        *g_x_Ptr = rotate[0][0]*origin_x + rotate[0][1]*origin_y + rotate[0][2]*origin_z;
	    *g_y_Ptr = rotate[1][0]*origin_x + rotate[1][1]*origin_y + rotate[1][2]*origin_z;
	    *g_z_Ptr = rotate[2][0]*origin_x + rotate[2][1]*origin_y + rotate[2][2]*origin_z;

        /* Update "gGammaIndex". */
        gGamma->gGammaIndex[index] = (*runID);
    }

    /* Copy gamma value to "result". */
    result[0] = *g_x_Ptr;
    result[1] = *g_y_Ptr;
    result[2] = *g_z_Ptr;
}

/* ================================================================ */
/* Parent function that calls "getGamma" to compute gamma value
   associated with object 1. */
void getGamma1(double *Gamma, double rotate[][3],int index, double result[])
{
    getGamma(&gGamma1, Gamma, rotate, index, result);
}


/* ================================================================ */
/* Parent function that calls "getGamma" to compute gamma value
   associated with object 2. */
void getGamma2(double *Gamma, double rotate[][3],int index, double result[])
{
    getGamma(&gGamma2, Gamma, rotate, index, result);
}


/* ================================================================ */
/*                                                             T
   This function generates returns the computed value of (gamma .P).
   To speed up the operation, the computed value is stored in the
   "gGammaP" array so that it can be referenced in future.  Further
   request for this computed value results in a simple retrieval
   operation instead of recomputation.


   Input Parameters
   ================
                                                    T
   *gGammaP     -- Structure to hold computed (gamma .p) information.
   *normals     -- Matrix containing the set of normals that defines the
                   faces that make up the object.
   *xBars       -- Matrix containing the coordinates of points on the
                   faces of the object corresponding to the normals
                   defined in the matrix "normals".
   rotate[3][3] -- Contains the transpose of the rotation matrix.
   p[3]         -- Positional vector of seed point.
   index        -- Index into "normals" and "xBars" matrices for which the
                                       T
                   corresponding (gamma .p) value is required.


   Output Parameters
   =================
                                   T
   result       -- Requested (gamma .p) value.
*/
__inline double getGammaP(struct gammaPData *gGammaP, struct gammaData *gGamma,
						  double *Gamma,
                          double rotate[][3], double p[], int index)
{
    /* Variable to store intermediate result. */
    double value[3];

    /* Pointer to storage. */
    double *gPtr;


    /* ====================== */
    /* Set "gGammaP" pointers. */
    gPtr = gGammaP->gGammaP + index;

    /*                                                    T
       Check to see if it is necessary to recompute (gamma .p).
    */
    if (gGammaP->gGammaPIndex[index] != (*runID))
    {
        /*                       T
           Need to compute (gamma .p).  Get gamma value first.
        */
        getGamma(gGamma, Gamma, rotate, index, value);

        /*               T
           Compute (gamma .p).
        */
        *gPtr = value[0] * p[0] + value[1] * p[1] + value[2] * p[2];

        /* Update "gGammaPIndex". */
        gGammaP->gGammaPIndex[index] = (*runID);
    }

    /*                        T
       Return result of (gamma .p).
    */
    return *gPtr;
}


/* ================================================================ */
/*                                                          T
   Parent function that calls "getGammaP" to compute (gamma1 .p1)
   value associated with object 1. */
double getGamma1P1(double *Gamma, double rotate[][3],double p[], int index)
{
    return getGammaP(&gGamma1P1, &gGamma1,Gamma, rotate, p, index);
}


/* ================================================================ */
/*                                                          T
   Parent function that calls "getGammaP" to compute (gamma1 .p2)
   value associated with object 1. */
double getGamma1P2(double *Gamma, double rotate[][3], double p[], int index)
{
    return getGammaP(&gGamma1P2, &gGamma1, Gamma, rotate, p, index);
}


/* ================================================================ */
/*                                                          T
   Parent function that calls "getGammaP" to compute (gamma2 .p1)
   value associated with object 2. */
double getGamma2P1(double *Gamma, double rotate[][3],
                   double p[], int index)
{
    return getGammaP(&gGamma2P1, &gGamma2, Gamma, rotate, p,index);
}


/* ================================================================ */
/*                                                          T
   Parent function that calls "getGammaP" to compute (gamma2 .p2)
   value associated with object 2. */
double getGamma2P2(double *Gamma, double rotate[][3],
                   double p[], int index)
{
    return getGammaP(&gGamma2P2, &gGamma2, Gamma, rotate, p,  index);
}


/* ================================================================ */
/* This function initializes the "gGammaIndex" array.

   Input Parameters
   ================
   *gGamma      -- Structure that holds computed gamma information.
*/
__inline void initgGammaIndex(struct gammaData *gGamma)
{
    int *i;

    for (i = gGamma->gGammaIndex; i < (gGamma->gGammaIndex + gGamma->size); i++)
        /* Zero element. */
        *i = INIT_ID;
}



/* ================================================================ */
/* This function initializes the "gGammaPIndex" array.

   Input Parameters
   ================
                                                       T
   *gGammaP     -- Structure that holds computed (gamma .p) information.
*/
__inline void initgGammaPIndex(struct gammaPData *gGammaP)
{
    int *i;

    for (i = gGammaP->gGammaPIndex; i < (gGammaP->gGammaPIndex + gGammaP->size); i++)
        /* Zero element. */
        *i = INIT_ID;
}



/* ================================================================ */
/* This function initializes the "gGamma?" data structure. */
void initGammaArrays(int size1, int size2, double *dWorkArray,
                     int *iWorkArray)
{
    runID = iWorkArray;

    gGamma1.gGamma = dWorkArray;
    gGamma1.gGammaIndex = runID + 1;
    gGamma1.size = size1;

    gGamma2.gGamma = gGamma1.gGamma + (size1 * 3);
    gGamma2.gGammaIndex = gGamma1.gGammaIndex + size1;
    gGamma2.size = size2;

    gGamma1P1.gGammaP = gGamma2.gGamma + (size2 * 3);
    gGamma1P1.gGammaPIndex = gGamma2.gGammaIndex + size2;
    gGamma1P1.size = size1;

    gGamma1P2.gGammaP = gGamma1P1.gGammaP + (size1 * 3);
    gGamma1P2.gGammaPIndex = gGamma1P1.gGammaPIndex + size1;
    gGamma1P2.size = size1;

    gGamma2P1.gGammaP = gGamma1P2.gGammaP + (size1 * 3);
    gGamma2P1.gGammaPIndex = gGamma1P2.gGammaPIndex + size1;
    gGamma2P1.size = size2;

    gGamma2P2.gGammaP = gGamma2P1.gGammaP + (size2 * 3);
    gGamma2P2.gGammaPIndex = gGamma2P1.gGammaPIndex + size2;
    gGamma2P2.size = size2;

    if ((*runID) == INT_MAX)
    {
        /* Initialize "runID". */
        (*runID) = FIRST_RUN_ID;

        /* Initialize data structure. */
        initgGammaIndex(&gGamma1);
        initgGammaIndex(&gGamma2);
        initgGammaPIndex(&gGamma1P1);
        initgGammaPIndex(&gGamma1P2);
        initgGammaPIndex(&gGamma2P1);
        initgGammaPIndex(&gGamma2P2);
    }
    else
        (*runID)++;
}
