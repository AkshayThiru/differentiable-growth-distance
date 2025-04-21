#ifndef __GAMMA_H
#define __GAMMA_H



/* Prototype declarations. */
void getGamma1(double *Gamma, double rotate[][3],int index, double result[]);

void getGamma2(double *Gamma, double rotate[][3],int index, double result[]);

double getGamma1P1(double *Gamma, double rotate[][3],double p[], int index);

double getGamma1P2(double *_Gamma, double rotate[][3], double p[], int index);

double getGamma2P1(double *Gamma, double rotate[][3], double p[], int index);

double getGamma2P2(double *Gamma, double rotate[][3], double p[], int index);


void initGammaArrays(int size1, int size2, double *dWorkArray,
                     int *iWorkArray);


#endif

