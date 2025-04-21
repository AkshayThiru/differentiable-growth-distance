#ifndef __SOLVER4c_H
#define __SOLVER4c_H


/* String needed by main driver routine to identify version of solver. */
#define SOLVER  "Solver 4c"


/* Prototype declarations. */
int solver(double *Gamma1, double *Gamma2,
		   int *adjSurf1, int *adjSurf2,
           int *adjFaceCount1, int *adjFaceCount2, int adjCol1,
           int adjCol2, double r1[][3], double r2[][3], double p1[],
           double p2[], int iIndex1[], int iIndex2[], int *nA, int *nB,
           int jIndex[][2], double result[], double dInv[][4], int dInvReq,int increment);
void get_origin_vector(double x[], double origin_x[], double r[][3]);


#endif

