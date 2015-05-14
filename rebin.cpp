#include "rebin.hpp"

/* to rebin input time-series */
/* consider input array as not that large */
off_t rebinning (float *x, off_t size, int nbins) {
 float sum;
 off_t newsize = size / nbins;

 for (off_t i=0, j, k=0; i<size; i+=nbins, k++) {
   sum = 0.;
   for (j=i; j<(i+nbins >= size ? size : i+nbins); j++) sum += x[j];
   x[k] = sum / nbins;
 }
 // zero-pad the remaining part of the input array
 for (off_t i=newsize; i<size; ++i) x[i] = 0.;

 return newsize;
}

/* to rebin input time-series */
/* using other array for output rebinned samples */
void rebinning_large (float *inp, float *x, unsigned long size, int nbins, off_t pos) {
 float sum;
 unsigned long newsize = size / nbins;

 for (unsigned long i=0, j, k=0; i<size; i+=nbins, k++) {
   sum = 0.;
   for (j=i; j<(i+nbins >= size ? size : i+nbins); j++) sum += inp[j];
   x[pos+k] = sum / nbins;
 } //for
}
