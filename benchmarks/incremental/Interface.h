#ifndef __INTERFACE_H
#define __INTERFACE_H



/* Prototype declaration. */
int gd(double *Gamma1, double *Gamma2,int *adjSurf1, int *adjSurf2,
       int *adjFaceCount1, int *adjFaceCount2, int adjCol1,
       int adjCol2, double r1[][3], double r2[][3], double p1[],
       double p2[], int jIndex[][2], double result[],
	   double dInv[][4], int dInvReq, int increment, int faces1,
	   int faces2, double *dWorkArray, int *iWorkArray);

void setup_gamma(double *normals, double *xBars, int sNum, double **Gamma);


#endif
