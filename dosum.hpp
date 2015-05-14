#ifndef _DOSUM_HPP
#define _DOSUM_HPP

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// 2-D arrays for FFA
extern float **y;
extern float **z;

/* implementation of the FFA algorithm */
void dosum (unsigned long P0, int k);

/*********************************************************************
 * starting filling of the y-array from intensity x-array (1st stage)
 * we have a 2-D array of M equations for P0 different arrival times
 **********************************************************************/
void array_filling (float *x, unsigned long P0, unsigned long M);

// with rebinning of folded profiles
void array_filling_rebin (float *x, unsigned long P0, unsigned long M, int mbins);

#endif //#ifndef _DOSUM_HPP
