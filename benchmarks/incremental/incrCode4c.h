#ifndef __INCRCODE4c_H
#define __INCRCODE4c_H

/* Prototype declarations. */
int incrSolver(double *Gamma1,double *Gamma2,
			   int *adjSurf1, int *adjSurf2,
               int *adjFaceCount1, int *adjFaceCount2,
               int adjCol1, int adjCol2, double r1[][3],
               double r2[][3], double p1[], double p2[],
               int iIndex1[], int iIndex2[], int *nA, int *nB,
               int jIndex[][2], double result[]);

#endif

