#ifndef _DOSNR_HPP
#define _DOSNR_HPP

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/******************************************************************************
* Function to calculate the maximum SNR
* for every M period from P0 to P0+1
* and maximum SNR for all of these periods
* Input args:
*  x - 2-D array with folded profiles
*  M - number of periods
*  P0 - base period
*  tres - time resolution
*  period_table - table with period, snr, number of snrs in that phase
*  index - current index of the table
******************************************************************************/
void dosnr (float **res, float **x, unsigned long M, unsigned long P0, int mbins, int pbins, double tres, double **snr, unsigned long **phase, double **period_table, unsigned long index, double snrmax_p[3][5], float **profiles);

/* sorting snr candidates in descending order */
void sort_snr (double **snr, unsigned long **phase, unsigned long M, unsigned long halfwidth);

/* sorting snrmax_p array */
void sort_snrmax (double snrmax_p[3][5], float **profiles, int size, int cols);

/* sorting period candidates in descending order */
void sort (double **x, unsigned long cols, unsigned long size, int sort_index);

/* sifting period candidates */
int sifting (double **x, double **tsift, unsigned long cols, unsigned long size, int sift_index, double tres, int mbins);

/* to rebin input time-series */
/* consider input array as not that large */
off_t extra_rebinning (float **x, unsigned long index, off_t size, int nbins);

/* to rebin input time-series */
/* using other array for output rebinned samples */
off_t extra_rebinning_large (float **res, float **x, unsigned long index, off_t size, int nbins);

/* to smooth the profiles (extra) */
/* in the end of array we will use the points from the beginning (wrapping) 
 *  * so the size will be the same */
off_t extra_smoothing_large (float **res, float **x, unsigned long index, off_t size, int nbins);

// the ratio of the pulse width to the period
// determined in ffasearch.hpp
extern double wpratio;

#endif //#ifndef _DOSNR_HPP
