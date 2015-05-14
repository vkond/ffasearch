/******************************************************************************
 * Function to calculate the maximum SNR
 * for every M period from P0 to P0+1
 * and maximum SNR for all of these periods
 * Input args:
 *  res - 2-D array with folded profiles
 *  x - 2-D auxiliary array with folded REBINNED profiles
 *  M - number of periods
 *  P0 - base period
 *  tres - time resolution
 *  period_table - table with period, snr, number of snrs in that phase
 *  index - current index of the table
 ******************************************************************************/
#include "dosnr.hpp"

void dosnr (float **res, float **x, unsigned long M, unsigned long P0, int mbins, int pbins, double tres, double **snr, unsigned long **phase, double **period_table, unsigned long index, double snrmax_p[3][5], float **profiles) {
 double aver, aver2, sigma;
 double sum, sum2;

 // halfwidth of effective pulse width to sift through the candidates
 unsigned long halfwidth = (unsigned long)ceil(P0 * wpratio * 0.5);

 unsigned long plength = P0 / mbins;
               halfwidth /= mbins;

 // extra rebinning
 if (pbins > 1) {
  for (unsigned long i=0; i<M; ++i)
//   extra_rebinning_large (res, x, i, plength, pbins);
   extra_smoothing_large (res, x, i, plength, pbins);

//  plength /= pbins;
//  halfwidth /= pbins;
 } // if (pbins > 1)

 // number of this phase in all M periods
 // but only for 3 largest snr
 unsigned long nphase[3] = { 0, 0, 0 };

 // clearing snr array
 for (long i=0; i<M; ++i) snr[0][i] = 0.;

 //first we will look for maximum and then calculate the aver and sigma in all array
 //except for 20% of points around the maximum (10% for every side)
 unsigned long shift = (unsigned long)ceil(0.1 * plength);

 // searching for maximum in the period
 for (unsigned long j=0; j<plength; ++j) for (unsigned long i=0; i<M; ++i) {
   phase[0][i] = i;
   if (x[j][i] > snr[0][i]) { snr[0][i] = x[j][i]; phase[1][i] = j; }
  } // for i
 for (unsigned long i=0; i<M; ++i) phase[2][i] = 1; // set "1" to all the fields

 // calculating aver, sigma and snr
 for (unsigned long i=0; i<M; ++i) {
  sum = 0.; sum2 = 0.;
  if (phase[1][i] > shift) for (unsigned long j=0; j<phase[1][i]-shift; ++j) { sum += x[j][i]; sum2 += x[j][i] * x[j][i]; }
  if (phase[1][i]+shift+1 < plength) for (unsigned long j=phase[1][i]+shift+1; j<plength; ++j) { sum += x[j][i]; sum2 += x[j][i] * x[j][i]; }

  aver   = sum  / ((phase[1][i] - shift > 0 ? phase[1][i] -shift : 0) + (phase[1][i] + shift + 1 < plength ? plength - phase[1][i] - shift - 1 : 0));
  aver2  = sum2 / ((phase[1][i] - shift > 0 ? phase[1][i] -shift : 0) + (phase[1][i] + shift + 1 < plength ? plength - phase[1][i] - shift - 1 : 0));
  sigma  = sqrt (aver2 - aver * aver);
  snr[0][i] = (snr[0][i] - aver) / sigma;
  snr[1][i] = aver;
  snr[2][i] = sigma;
 } // for i

  // looking for 3 best candidates
  sort_snr (snr, phase, M, halfwidth);
  // computing nphase
  for (unsigned long j=0; j<M; ++j) if (phase[1][j] == phase[1][0]) nphase[0]++;
  if (phase[1][1] == phase[1][0]) nphase[1] = nphase[0]; else for (unsigned long j=0; j<M; ++j) if (phase[1][j] == phase[1][1]) nphase[1]++;
  if (phase[1][2] == phase[1][1]) nphase[2] = nphase[1]; 
   else { if (phase[1][2] == phase[1][0]) nphase[2] = nphase[0]; 
           else for (unsigned long j=0; j<M; ++j) if (phase[1][j] == phase[1][2]) nphase[2]++;
   }

  bool is_exist = false; // if candidate is already in the list
  long cand_index = -1; // index of the candidate which is already in the list
  // filling 3 best candidates
  for (int y=0; y<3; y++) {
   period_table[0][index+y] = ((double)P0 + ((double)phase[0][y] / (double)(M-1))) * tres * 0.001; // period of the candidate in ms
   period_table[1][index+y] = snr[0][y];
   period_table[2][index+y] = ((double)phase[1][y] / (double)plength) * 360.; // phase in degrees 
   period_table[3][index+y] = (double)nphase[y];
   period_table[4][index+y] = (double)pbins;

   is_exist = false;
   // checking if the candidate is already in the list
   for (int u=0; u<3; ++u) {
    if (phase[1][y] >= (long)snrmax_p[u][2]-halfwidth && phase[1][y] <= (long)snrmax_p[u][2]+halfwidth) { 
     is_exist = true;
     cand_index = u;
     break;
    } // if
   } // for u
   if (is_exist) { // if this candidate exists already
    if (snr[0][y] > snrmax_p[cand_index][0]) {
     snrmax_p[cand_index][0] = snr[0][y];
     snrmax_p[cand_index][1] = plength;
     snrmax_p[cand_index][2] = ((double)phase[1][y] / (double)plength) * 360.; // phase in degrees
     snrmax_p[cand_index][3] = ((double)P0 + ((double)phase[0][y] / (double)(M-1))) * tres * 0.001; // period of the candidate in ms
     snrmax_p[cand_index][4] = (double)pbins;
     for (long j=0; j<plength; ++j) profiles[cand_index][j] = (float)(((double)x[j][phase[0][y]] - snr[1][y]) / snr[2][y]);
    }
   } else {
    for (int u=0; u<3; ++u) { // forming three best candidates' profiles
     if (snr[0][y] > snrmax_p[u][0]) {
      snrmax_p[2][0] = snr[0][y];
      snrmax_p[2][1] = plength;
      snrmax_p[2][2] = ((double)phase[1][y] / (double)plength) * 360.; // phase in degrees
      snrmax_p[2][3] = ((double)P0 + ((double)phase[0][y] / (double)(M-1))) * tres * 0.001; // period of the candidate in ms
      snrmax_p[2][4] = (double)pbins;
      for (long j=0; j<plength; ++j) profiles[2][j] = (float)(((double)x[j][phase[0][y]] - snr[1][y]) / snr[2][y]);
      break;
     }
    } // for u
   } // else if is_exist
   sort_snrmax (snrmax_p, profiles, 3, 5);
  } // for y
}

/* sorting snr candidates in descending order just for 3 best candidates */
void sort_snr (double **snr, unsigned long **phase, unsigned long M, unsigned long halfwidth) {
 double max, tempd;
 unsigned long jmax, templ;
 for (unsigned long i=0; i<(M-1 >= 3 ? 3 : M-1); ++i) {
  max = snr[0][i];
  jmax = i;
  for (unsigned long j=i; j<M; ++j) { //look for the maximum
   if (snr[0][j] > max) { if (phase[2][j] == 1) { max = snr[0][j]; jmax = j; } }
  } //for j
  // swap max element and current i-element
  for (unsigned long k=0; k<3; ++k) { 
   tempd = snr[k][i]; snr[k][i] = snr[k][jmax]; snr[k][jmax] = tempd;
   templ = phase[k][i]; phase[k][i] = phase[k][jmax]; phase[k][jmax] = templ; 
  } // for k
  // marking the phases that falls into the halfwidth from the maximum
  // !!!!!!!!!!!!!! why I am doing this?? I don't use it at all... I guess, I have changes my mind at some point
  // It seems like I don't need this phase[2][]
  for (unsigned long k=i+1; k<M; ++k) {
   if (phase[1][k] >= phase[1][i]-halfwidth && phase[1][k] <= phase[1][i]+halfwidth) phase[2][k] = 0;
  }
 } //for i
}

/* sorting period candidates in descending order */
void sort (double **x, unsigned long cols, unsigned long size, int sort_index) {
 double max, temp;
 unsigned long jmax;
 for (unsigned long i=0; i<size-1; ++i) {
  max = x[sort_index][i];
  jmax = i;
  for (unsigned long j=i; j<size; ++j) { //look for the maximum
   if (x[sort_index][j] > max) { max = x[sort_index][j]; jmax = j; }
  } //for j
  // swap max element and current i-element
  for (unsigned long k=0; k<cols; ++k) { temp = x[k][i]; x[k][i] = x[k][jmax]; x[k][jmax] = temp; }
 } //for i
}

/* sifting period candidates */
int sifting (double **x, double **tsift, unsigned long cols, unsigned long size, int sift_index, double tres, int mbins) {

 // halfwidth of effective pulse width to sift through the candidates
// unsigned long halfwidth;
 double halfwidth = 0.5 * wpratio * 360.; // in degrees

 int sift_size = 0;
 // we will sift also looking for _exact_ periods which are possible due to the range of different pbins but having different SNRs. 
 // also if phases of candidates are close than we will leave the first (i.e. with highest SNRs)
 // first candidates are always with higher SNRs because there was a sorting before
 bool phase_exist = false; // flag - if phase is already exist
 // main loop
 for (unsigned long i=0; i<size; ++i) {
  phase_exist = false;
//  halfwidth = (unsigned long)ceil(floor((x[0][i] * 1000.)/(tres * mbins * x[4][i])) * wpratio * 0.5);
  for (unsigned long j=0; j<sift_size; ++j) {
   if (x[0][i] == tsift[0][j]) { phase_exist = true; tsift[3][j] += x[3][i]; break; } // periods are the same
//   if ((long)x[sift_index][i] >= (long)tsift[sift_index][j]-halfwidth && (long)x[sift_index][i] <= (long)tsift[sift_index][j]+halfwidth) 
   if (x[sift_index][i] >= tsift[sift_index][j]-halfwidth && x[sift_index][i] <= tsift[sift_index][j]+halfwidth) 
    { phase_exist = true; tsift[3][j] += x[3][i]; break; }
  } // for j
  if (!phase_exist) {
   for (int m=0; m<cols; m++) tsift[m][sift_size] = x[m][i];
   sift_size++;
  } 
 } //for i

 return sift_size;
}


/* sorting snrmax_p array */
void sort_snrmax (double snrmax_p[3][5], float **profiles, int size, int cols) {
 double max, temp;
 int jmax;
 float temp1;
 long pmax;
 for (int i=0; i<size-1; ++i) {
  max = snrmax_p[i][0];
  jmax = i;
  for (int j=i; j<size; ++j) { // look for the maximum
   if (snrmax_p[j][0] > max) { max = snrmax_p[j][0]; jmax = j; }
  } // for j
  // swap max element and current i-element
  for (int k=0; k<cols; ++k) { temp = snrmax_p[i][k]; snrmax_p[i][k] = snrmax_p[jmax][k]; snrmax_p[jmax][k] = temp; }
  pmax = (unsigned long)snrmax_p[i][1];
  if ((unsigned long)snrmax_p[jmax][1] > pmax) pmax = (unsigned long)snrmax_p[jmax][1];
  for (unsigned long t=0; t<pmax; ++t) { temp1 = profiles[i][t]; profiles[i][t] = profiles[jmax][t]; profiles[jmax][t] = temp1; }
 } // for i
}

/* to rebin input time-series */
/* consider input array as not that large */
off_t extra_rebinning (float **x, unsigned long index, off_t size, int nbins) {
 float sum;
 off_t newsize = size / nbins;

 for (off_t i=0, j, k=0; i<size; i+=nbins, k++) {
  sum = 0.;
  for (j=i; j<(i+nbins >= size ? size : i+nbins); j++) sum += x[j][index];
  x[k][index] = sum / nbins;
 }

 return newsize;
}

/* to rebin input time-series */
/* using other array for output rebinned samples */
off_t extra_rebinning_large (float **res, float **x, unsigned long index, off_t size, int nbins) {
 float sum;
 unsigned long newsize = size / nbins;

 for (unsigned long i=0, j, k=0; i<size; i+=nbins, k++) {
  sum = 0.;
  for (j=i; j<(i+nbins >= size ? size : i+nbins); j++) sum += res[j][index];
  x[k][index] = sum / nbins;
 } //for

 return newsize;
}

/* to smooth the profiles (extra) */
/* in the end of array we will use the points from the beginning (wrapping) 
 * so the size will be the same */
off_t extra_smoothing_large (float **res, float **x, unsigned long index, off_t size, int nbins) {
 float sum = 0., tmp;

 for (off_t r=0; r<nbins; ++r) sum += res[r][index];
 tmp = res[0][index];
 x[0][index] = sum / nbins;

 for (off_t t=0; t<size - nbins; ++t) {
  sum += (res[nbins+t][index] - tmp);
  tmp = res[t+1][index];
  x[t+1][index] = sum / nbins;
 }

 for (off_t s=size-nbins, i=0; s<size-1; ++s, ++i) {
  sum += (res[i][index] - tmp);
  tmp = res[s+1][index];
  x[s+1][index] = sum / nbins;
 }
}
