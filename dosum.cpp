/****************************************************************************
 * This file holds the implementation of the base function of FFA algorithm
 * to make efficient summing of the input array for the range of periods
 * Input arguments are:
 * - array
 * - base period P0
 * - k, power of 2 which is 2^k = M - number of periods
 * Output is 2-D array which are folded M profiles with
 * base period P0 (maybe I will change this if this array gets very huge and
 * then I will call for function calculating the SNR)
 * size/P0 = M - should be
 * Here I suggest, that size is small enough for the array to be allocated at once
 *  otherwise, I will need to rewrite this function in different manner
******************************************************************************/

#include "dosum.hpp"

void dosum (unsigned long P0, int k) {
 
 unsigned long M = 1 << k; // number of periods

 // temporary pointers
 float **ptrz, **ptry;
 unsigned long shift, shift1, st, s, m;
 // transforming the 2-D array of sums (next stages)
 for (st=1; st<k; ++st) { // this loop is loop on stages, k is 2^k = M
  if (st%2 == 0) { ptrz = y; ptry = z; } else { ptrz = z; ptry = y; }
  shift = 1 << st;
  shift1 = 2 * shift;
  for (m=0; m<P0; ++m) {       // loop on every sample (count) in period
   for (s=0; s<M; s+=shift1) { // loop on the groups of similar additions 
                               // for example, on the first stage the index s runs for every sum B0

    // these two loops are on sums withis one group
    // for the first stage these loops are for B0 and B1, B2 and B3
    // for the second stage, for C0 and C1, C2 and C3, C4 and C5 and so on
    for (unsigned long n=s, j=m, w=s, v=s+shift; n<s+shift1; n+=2, ++j, ++w, ++v) { if (j==P0) j-=P0;
     ptry[m][n] = ptrz[m][w] + ptrz[j][v];
    }
    for (unsigned long n=s+1, j=m+1, w=s, v=s+shift; n<s+shift1; n+=2, ++j, ++w, ++v) { if (j==P0) j-=P0;
     ptry[m][n] = ptrz[m][w] + ptrz[j][v];
    }
   } // for s
  } //for m
 } // for st
}

/*********************************************************************
 * starting filling of the y-array from intensity x-array (1st stage)
 * we have a 2-D array of M equations for P0 different arrival times
**********************************************************************/

void array_filling (float *x, unsigned long P0, unsigned long M) {
  unsigned long p1 = P0-1, p2 = 2*P0;
  for (unsigned long m=0; m<p1; ++m) for (unsigned long n=0, np=m, np1=m+P0; n<M; n+=2, np+=p2, np1+=p2) {
     z[m][n]   = x[np] + x[np1];
     z[m][n+1] = x[np] + x[np1+1];
   } // for n
   for (unsigned long n=0, np=p1, np1=P0; n<M; n+=2, np+=p2, np1+=p2) {
     z[p1][n]   = x[np] + x[np1+p1];
     z[p1][n+1] = x[np] + x[np1];
   }
}

void array_filling_rebin (float *x, unsigned long P0, unsigned long M, int mbins) {
  unsigned long np, np1, plength = P0/mbins;
  float x0, x1, x2;
  for (unsigned long m=0; m<plength-1; ++m) for (unsigned long n=0; n<M; n+=2) {
     np = n * P0 + m * mbins;
     np1 = (n+1) * P0 + m * mbins;
     x0 = 0.; x1 = 0.; x2 = 0.;
     for (int l=np, o=np1, p=np1+mbins; l<np+mbins; ++l, ++o, ++p) { x0 += x[l]; x1 += x[o]; x2 += x[p]; }
     x0 /= mbins; x1 /= mbins; x2 /= mbins;
     z[m][n]   = x0 + x1;
     z[m][n+1] = x0 + x2;
   } // for n
   for (unsigned long n=0; n<M; n+=2) {
     np = n * P0 + P0 - mbins;
     np1 = (n+1) * P0 + P0 - mbins;
     x0 = 0.; x1 = 0.; x2 = 0.;
     for (int l=np, o=np1, p=np1+mbins-P0; l<np+mbins; ++l, ++o, ++p) { x0 += x[l]; x1 += x[o]; x2 += x[p]; }
     x0 /= mbins; x1 /= mbins; x2 /= mbins;
     z[plength-1][n]   = x0 + x1;
     z[plength-1][n+1] = x0 + x2;
   } // for n
}
