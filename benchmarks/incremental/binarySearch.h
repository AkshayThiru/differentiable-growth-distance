#ifndef __BINARY_SEARCH_H
#define __BINARY_SEARCH_H



/* Prototype declarations. */
int binarySearch(void (*getGamma)(), double *Gamma,
                 double rotate[][3], double deltaX[],
                 double xMinusP[], double sigma, int *adjSurf,
                 int adjFaceCount, int *activeFace,
                 double *stepSize);
int binarySearchFinal(void (*getGamma)(), double *Gamma,
                      double rotate[][3], double deltaX[], 
                      double xMinusP[], double sigma, double descent,
                      int *adjSurf, int adjFaceCount, int *activeFace,
                      double *stepSize);


#endif

