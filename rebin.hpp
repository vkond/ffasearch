#ifndef _REBIN_HPP
#define _REBIN_HPP

#include <stdlib.h>

/* to rebin input time-series */
/* consider input array as not that large */
off_t rebinning (float *x, off_t size, int nbins);

/* to rebin input time-series */
/* using other array for output rebinned samples */
void rebinning_large (float *inp, float *x, unsigned long size, int nbins, off_t pos);

#endif //#ifndef _REBIN_HPP
