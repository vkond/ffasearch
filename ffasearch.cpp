/*****************************************************************
 *
 *          Search time series with FFA algorithm
 *
 * written by Vlad Kondratiev, Dec 2006, West Virginia University
 *
 *****************************************************************/
#include "ffasearch.hpp"

// get time in readable format, days, hours, mins and seconds
char *gettime (time_t seconds, char *dest) {
 int days = (int)(seconds / 86400);
 int hours = (int)((seconds - days * 86400) / 3600);
 int minutes = (int)((seconds - days * 86400 - hours * 3600) / 60);
 int secs = (int)(seconds - days * 86400 - hours * 3600 - minutes * 60);
 sprintf (dest, "%d days %d hours %d minutes %d seconds", days, hours, minutes, secs);
 return dest;
}

// get closest power of 2
int get_power2 (off_t nsamp, long p0) {
 off_t size = nsamp / (off_t)p0 + 1;
 int n = 0;
 while ((1 << (++n)) < size) ; 
 return --n;
}

// number of periods to process
unsigned long number_of_periods (unsigned long p1, unsigned long p2, off_t Nsamples, int mbins, bool is_period_res) {
 unsigned long Nper = 0;
 if (!is_period_res) {
  int k;
  unsigned long M;
  for (unsigned long p=p1; p<=p2; p++) {
    k = get_power2 (Nsamples, p * mbins);
    if (is_zero_padding) k += 1; // use zero-padding
    M = 1 << k;
    Nper += M;
   } //for p
 } else { // is_period_res = true
   for (unsigned long p=p1; p<=p2; p++) Nper += Nsamples;
 }
 return Nper;
}

// calculating the maximum value of M and p0 to allocate at once arrays of largest necessary size
// to save time for reallocating the arrays
// also calculating the usual pulse width as 3% of the lowest period _plow_
// and suggesting the use of rebinning if is needed
void planning (double plow, double phigh, double tres, int nbins, int mbins, int pbins, off_t Nsamples, bool is_period_res, bool is_zero_padding, unsigned long& Mmax, unsigned long& Mp0max) {

 unsigned long p0low = (unsigned long)floor((plow * 1000.)/(tres * nbins * mbins));
 unsigned long p0high = (unsigned long)floor((phigh * 1000.)/(tres * nbins * mbins));
 unsigned long p0; // base period

 if (!is_period_res) {
  int k;            // power of 2, 2^k = M
  unsigned long M;  // number of base periods P0 in time series
  // loop on increment of p0
  for (p0=p0low; p0<=p0high; p0++) {
   k = get_power2 (Nsamples/nbins, p0 * mbins);
   if (is_zero_padding) k += 1; // use zero-padding
   M = 1 << k;
   if (M < 2) { // FFA can not work if M < 2 because sums shifted at least on time in 1 samples are required
    printf ("# - warning -\n");
    printf ("# value of largest period is too high for used Nsamples!\n");
    printf ("# decrease _phigh_ or increase Nsamples!\n");
    printf ("# FFA can not work if M < 2!\n");
    exit (1);
   }
   if (p0 == p0low) { Mmax = M; Mp0max = M * p0 * mbins; }
   if (M > Mmax) Mmax = M;
   if (M * p0 > Mp0max) Mp0max = M * p0 * mbins;
  } // for p0
  if ((Nsamples/nbins) > Mp0max) Mp0max = (unsigned long)(Nsamples/nbins);
 } else { Mmax = Nsamples; Mp0max = Mmax * p0high * mbins; }  // when is_period_res = true, Nsamples = M is the planning function parameters
 
 double effwidth = wpratio * plow;
 unsigned long factor = (unsigned long)rint((effwidth * 1000.) / tres); 
 printf ("# probable lowest effective pulse width = %.1lf ms\n", effwidth);
 printf ("# suggested rebin factor = %ld\n", factor);
 printf ("#    raw-data decimated by %d, folded profiles rebinned by %d, extra rebinning by %d  -> \n", nbins, mbins, pbins);
 printf ("#    total rebin factor = %d\n", nbins * mbins * pbins);
 if (factor > nbins * mbins * pbins) {
  printf ("# - warning - \n");
  printf ("# you may want to increase your rebin factor for more efficiency!\n");
 }
 printf ("#    covered range of duty cycles is %.1f - %.1f\%\n", ((tres * nbins * mbins)/phigh)*0.1, ((tres * nbins * mbins * pbins)/plow)*0.1);

 if (((is_period_res ? Nsamples * p0high * mbins : Nsamples)/nbins)*sizeof(float) > block) {
  printf ("# - WARNING - \n");
  printf ("# Original number of samples is %ld with downsampling by %d\n", is_period_res ? (unsigned long)(Nsamples * p0high * mbins * nbins) : (unsigned long)Nsamples, nbins);
  printf ("# Number of samples %ld is too large!\n", (unsigned long)((is_period_res ? Nsamples * p0high * mbins : Nsamples)/nbins));
  printf ("# It is larger than maximum block size of %ld bytes\n", (unsigned long)block);
  printf ("# Use downsampling or smaller window instead\n");
  printf ("#\n");
  exit (1);
 }
}

/*********************** M  A  I  N ********************************/
int main (int argc, char *argv[]) {

 if (argc == 1) help (argv[0]);

 char *outname, *timepres, *prespres, *elapsed_time;
 outname = new char[255];
 output_dir = new char[255];
 filestem = new char[255];
 if (!outname || !output_dir || !filestem) { perror("new output_dir | outname | filestem"); exit (1); }
 timepres = new char[255];
 prespres = new char[255];
 elapsed_time = new char[255];
 if (!timepres || !prespres || !elapsed_time) { perror("new timepres | prespres | elapsed_time"); exit (1); }
 memset (outname, 0, 255);
 memset (timepres, 0, 255);
 memset (prespres, 0, 255);
 memset (elapsed_time, 0, 255);
 memset (output_dir, 0, 255);
 memset (filestem, 0, 255);
 sprintf (output_dir, "%s", ".");
 sprintf (filestem, "%s", "\0");
 
 int next = parse_command_line (argc, argv);
 if (argc < next+1) { printf ("Not point the files to process!\n"); exit (1); }

 struct stat info;
 if (stat(argv[next], &info) != 0) { perror("stat"); exit(1); }

 if (is_sigproc) { // if tim-file is in sigproc format
    sigobj.initialize (argv[next]);
    if (!sigobj.is_initialized()) {
     printf ("Bad sigproc header format in the input file %s!", argv[next]);
     exit (1);
    } else { // if header is good
      header = sigobj.get_header_size();
      tres = 1.e6 * sigobj.get_tsamp(); // converting them to mcs
      dm = sigobj.get_refdm();
    }
 }
 filesamples = (info.st_size - header) / sizeof(float);

 if (Nsamples == 0) Nsamples = filesamples;
 if (left_edge > filesamples-1) left_edge = 0;
 if (left_edge + Nsamples > filesamples) Nsamples = filesamples - left_edge;

 if (nbins >= Nsamples) { nbins = 1; } // if nbins is too large

 unsigned long M;  // number of base periods P0 in time series
 // if period_res was set in command line
 if (is_period_res) {
  // if mbins is set reduce period_res by mbins
  if (period_res/mbins > tres * nbins) is_period_res = false;
   else {
    M = 1 << (get_power2 ((off_t)((tres * nbins)/(period_res / mbins)), 1));
    if (M < 2) is_period_res = false;
     else { period_res = (tres * nbins / M) * mbins;
            Nsamples = (M * ((unsigned long)floor((phigh * 1000.)/(tres*nbins)))) * nbins;
            if (left_edge + Nsamples > filesamples || nbins >= Nsamples) is_period_res = false;
	  } 
   }
 }
 if (is_Mcmd) {
  Mcmd = 1 << (get_power2 (Mcmd, 1));
  if (Mcmd >= 2) { if ((is_period_res && Mcmd < M) || !is_period_res) { 
	     is_period_res = true;
	     M = Mcmd; 
	     period_res = (tres * nbins / M) * mbins;
	     Nsamples = (M * ((unsigned long)floor((phigh * 1000.)/(tres*nbins)))) * nbins; 
	     if (left_edge + Nsamples > filesamples || nbins >= Nsamples) is_period_res = false;
          }
   } // Mcmd >= 2
 } // if is_Mcmd

 // checking mbins
 // it should not be larger than half of the lowest possible period (in samples)
 if (mbins > ((unsigned long)floor((plow * 1000.)/(tres*nbins)))/2) mbins = 1;

 // doing the same for pbins with already applied mbins
 // it should not be larger than half of the lowest possible period (in samples)
 if (pbins > ((unsigned long)floor((plow * 1000.)/(tres*nbins*mbins)))/2) pbins = 1;

 if (!is_period_res) sprintf (prespres, "period resolution [%s]: %lg  -  %lg ms [to process %ld periods]", "auto", (tres * nbins * mbins * 0.001) / (1 << ((is_zero_padding ? 1 : 0) + get_power2 (Nsamples/nbins, (unsigned long)floor((plow * 1000.)/(tres*nbins))))), (tres * nbins * mbins * 0.001) / (1 << ((is_zero_padding ? 1 : 0) + get_power2 (Nsamples/nbins, (unsigned long)floor((phigh * 1000.)/(tres*nbins))))), number_of_periods ((unsigned long)floor((plow * 1000.)/(tres*nbins*mbins)), (unsigned long)floor((phigh * 1000.)/(tres*nbins*mbins)), Nsamples/nbins, mbins, false));
  else sprintf (prespres, "period resolution [%s]: %lg ms [to process %ld periods]", "fixed", period_res * 0.001, number_of_periods ((unsigned long)floor((plow * 1000.)/(tres*nbins*mbins)), (unsigned long)floor((phigh * 1000.)/(tres*nbins*mbins)), M, mbins, true));

 printf ("# File: %s\n", argv[next]);
 printf ("# Header = %ld bytes  Total samples = %ld  sigproc format - %s\n", header, (unsigned long)filesamples, (is_sigproc ? "yes" : "no"));
 printf ("# Left pos = %ld   Nsamples = %ld   Block size = %ld bytes\n", (unsigned long)left_edge, (unsigned long)Nsamples, (unsigned long)block);
 printf ("# Original time resolution = %lg mcs\n", tres);
 printf ("# DM = %lg pc/cm^3\n", dm);
 printf ("# Decimating by %ld   [used tres: %lg mcs] [used Nsamples: %ld]\n", nbins, nbins * tres, Nsamples/nbins);
 printf ("# Zero-padding the data... %s\n", is_zero_padding ? "yes" : "no");
 printf ("# Search in period range: %lg - %lg ms\n", plow, phigh);
 printf ("#      %s\n", prespres);
 printf ("#      rebinning folded profiles by %ld in the possible range [%ld - %ld]\n", mbins, 1, ((unsigned long)floor((plow * 1000.)/(tres*nbins)))/2);
 printf ("#      largest extra rebinning of folded profiles by %ld in the possible range [%ld - %ld]\n", pbins, 1, ((unsigned long)floor((plow * 1000.)/(tres*nbins*mbins)))/2);
 printf ("#\n");
 printf ("# output directory: %s\n", output_dir);
 printf ("# filestem:  %s    candidate file:  %s%s\n", (strcmp(filestem, "\0") == 0 ? "undefined" : filestem), (strcmp(filestem, "\0") == 0 ? basename(argv[next]) : filestem), ".ffa");
 printf ("# candidate file to be * %s *\n", is_to_append ? "APPENDED" : "OVERWRITTEN");
 printf ("#\n");
 if (block > MAXBLOCKSIZE) { 
  printf ("# - warning - \n"); 
  printf ("# block size of %ld bytes is too large!\n", (unsigned long)block);
  printf ("# hope, you know what are you doing...\n");
  printf ("#\n");
 }

 // calculating the maximum value of M and p0 to allocate at once arrays of largest necessary size
 // to save time for reallocating the arrays
 unsigned long Mmax, Mp0max;
 printf ("# planning ...\n");
 planning (plow, phigh, tres, nbins, mbins, pbins, is_period_res ? M : Nsamples, is_period_res, is_zero_padding, Mmax, Mp0max);
 printf ("# ...done\n");

 float *x;
 x = new float[Mp0max];
 if (!x) { perror("new x"); exit (1); }
 memset (x, 0, Mp0max * sizeof(float));

 printf ("# reading... "); fflush (NULL);
 unsigned long reading_step = ((unsigned long)((Nsamples*sizeof(float) > block ? block/sizeof(float) : Nsamples) / nbins)) * nbins; // we will read the file with this step in samples
 int in;
 off_t pos = header + left_edge * sizeof(float); // to read from this position in bytes
 off_t x_index = 0;  // position in x-array to write binned data
 int bytesread;      // number of bytes were read

 float *inp;
 inp = new float[reading_step];
 if (!inp) { perror("new inp"); exit (1); }

 while (x_index < Mp0max) {
  if ((in = open (argv[next], O_RDONLY | O_LARGEFILE)) == -1) { perror("open"); exit(1); };
  lseek (in, pos, SEEK_SET);
  memset (inp, 0, reading_step * sizeof(float));
  bytesread = read (in, inp, reading_step * sizeof(float));
  close (in);

  // downsampling
  if (bytesread == reading_step * sizeof(float)) {
   if (nbins > 1) rebinning_large (inp, x, reading_step, nbins, x_index);
    else for (off_t a=0; a<reading_step; ++a) x[x_index+a] = inp[a];
  }  

  x_index += (reading_step / nbins);
  pos += reading_step * sizeof(float);
 } // while pos

 delete (inp);
 printf ("done\n");
 printf ("#\n");

 // change time res also if downsampling is used
 // also change Nsamples 
 if (nbins > 1) { tres *= nbins; Nsamples /= nbins; }

 unsigned long p0low = (unsigned long)floor((plow * 1000.)/(tres * mbins));
 unsigned long p0high = (unsigned long)floor((phigh * 1000.)/(tres * mbins));
 unsigned long p0; // base period
 int k;            // power of 2, 2^k = M
 // maximum number of periods to process
 unsigned long Nperiods_max = number_of_periods (p0low, p0high, is_period_res ? M : Nsamples, mbins, is_period_res);

 // this array keeps SNR, aver, and sigma for every candidate of single FFA transaction (for specific P0)
 double **snr = (double **) calloc (3, sizeof(double *));
 if (!snr) { perror("calloc snr"); exit (1); }
 for (long i=0; i<3; ++i) {
  snr[i] = (double *) calloc (Mmax, sizeof (double));
  if (!snr[i]) { perror("calloc snr[i]"); exit (1); }
 }

 // phase[0][] keeps the value "i", which is p = p0 + i*(tres/(M-1)) to calculate exact period
 // phase[1][] keeps the phase in the folded profile
 // phase[2][] keeps 0 or 1. If candidate is not connected with previously detected one (i.e. it has the phase
 //   not to close with better candidate) than this cell has "1", otherwise "0"
 //   this is necessary to sift the candidates
 unsigned long **phase = (unsigned long **) calloc (3, sizeof(unsigned long));
 if (!phase) { perror("calloc phase"); exit (1); }
 for (int e=0; e<3; e++) {
  phase[e] = (unsigned long *) calloc (Mmax, sizeof (unsigned long));
  if (!phase[e]) { perror("calloc phase[e]"); exit (1); }
 }

 // 2-D array to hold info about:
 // 1st col - period
 // 2nd col - snr
 // 3rd col - phase
 // 4th col - nphase
 // 5th col - specific pbin
 double **period_table;
 unsigned long prange = p0high-p0low+1;
 period_table = new double*[5];
 if (!period_table) { perror("new period_table"); exit (1); }
 for (int h=0; h<5; ++h) {
  period_table[h] = new double[3*pbins*prange]; // 3 times larger because we want to keep 3 max snr from every M 
                                                // and then pbins larger because we are doing the same for the range of pbins (from 1 to pbins)
  if (!period_table[h]) { perror("new period_table[h]"); exit (1); }
  memset (period_table[h], 0, 3 * prange * pbins * sizeof(double));
 }

 // it is more efficient to have the arrays with size of power of 2
 unsigned long psize = 1 << (1 + get_power2(p0high, 1));

 // array that keeping 3 candidates with highest snrs during the processing
 // and their P0, and their phase, and their period, and their pbin: 
 // [0][0] - snr of first candidate, [0][1] - its P0, [0][2] - its phase, [0][3] - its period, [0][4] - its pbin
 double snrmax_p[3][5] = { {0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0.} };
 // array of 3 best profiles
 float **profiles = (float **) calloc (3, sizeof(float *));
 if (!profiles) { perror("calloc profiles"); exit (1); }
 for (long i=0; i<3; ++i) {
  profiles[i] = (float *) calloc (psize, sizeof (float));
  if (!profiles[i]) { perror("calloc profiles[i]"); exit (1); }
 }

 // 2-D arrays for doing FFA
 y = (float **)calloc (psize, sizeof(float *));
 if (!y) { perror("calloc y"); exit (1); }
 for (long i=0; i<psize; ++i) {
  y[i] = (float *)calloc (Mmax, sizeof (float));
  if (!y[i]) { perror("calloc y[i]"); exit (1); }
 }

 z = (float **)calloc (psize, sizeof(float *));
 if (!z) { perror("calloc z"); exit (1); }
 for (long i=0; i<psize; ++i) {
  z[i] = (float *)calloc (Mmax, sizeof (float));
  if (!z[i]) { perror("calloc z[i]"); exit (1); }
 }

 time_t starttime = time(NULL);
 fflush (NULL);
 if (is_verbose) { printf ("# processed: %3d%%", 0); fflush (NULL); }
 unsigned long index = 0, Np = 0;
 // loop on increment of p0
 for (p0=p0low; p0<=p0high; p0++, index++) {
  k = get_power2 (is_period_res ? M * p0 * mbins : Nsamples, p0 * mbins);
  if (!is_period_res) {
   if (is_zero_padding) k += 1; // use zero-padding
   M = 1 << k;
  }

  // starting filling of the z-array from intensity x-array (1st stage)
  if (mbins > 1) array_filling_rebin (x, p0 * mbins, M, mbins); 
   else array_filling (x, p0, M);
  // search
//  dosum (p0/mbins, k);
  dosum (p0, k);
  // calc snr
  if (k%2 == 0) { 
    for (int er = 1; er<=pbins; ++er) 
     // result is in y-array
     // z-array is auxiliary
     dosnr (y, z, M, p0 * mbins, mbins, er, tres, snr, phase, period_table, 3*pbins*index + 3*(er-1), snrmax_p, profiles);
//     dosnr (y, z, M, p0 * mbins, mbins, pbins, tres, snr, phase, period_table, 3*pbins*index, snrmax_p, profiles);
  } else {
    for (int er = 1; er<=pbins; ++er) 
     // result is in z-array
     // y-array is auxiliary
     dosnr (z, y, M, p0 * mbins, mbins, er, tres, snr, phase, period_table, 3*pbins*index + 3*(er-1), snrmax_p, profiles);
//     dosnr (z, y, M, p0 * mbins, mbins, pbins, tres, snr, phase, period_table, 3*pbins*index, snrmax_p, profiles);
  }

  if (is_verbose) { Np += number_of_periods(p0, p0, is_period_res ? M : Nsamples, mbins, is_period_res); printf ("\b\b\b\b%3d%%", (Np*100)/Nperiods_max); fflush (NULL); }

 } //for p0

 if (is_verbose) printf ("\n");
 time_t endtime = time(NULL);
 sprintf (elapsed_time, "Time elapsed: %ld seconds (%s)", (unsigned long)(endtime - starttime), gettime(endtime - starttime, timepres));
 printf ("# %s\n", elapsed_time);
 printf ("#\n");

 // forming the name of the output file
 sprintf (outname, "%s/%s%s", output_dir, (strcmp(filestem, "\0") == 0 ? basename(argv[next]) : filestem), ".ffa");

 // creating current postscript file
 plotting (period_table, 3*pbins*prange, outname, plow, phigh, prespres, elapsed_time, argc, argv, next, tres, nbins, mbins, pbins, dm, is_zero_padding, snrmax_p, profiles, sigobj);

 // sorting period_table
 sort (period_table, 5, 3 * pbins * prange, 1);

 // writing all period candidates into the file
 FILE *candlist;
 if ((candlist = fopen(outname, is_to_append ? "at" : "wt")) == NULL) { perror("fopen"); exit (1); }
 fprintf (candlist, "# File: %s\n", argv[next]);
 fprintf (candlist, "# DM = %lg pc/cm^3\n", dm);
 fprintf (candlist, "# Decimating by %ld   [used tres: %lg mcs] [used Nsamples: %ld]\n", nbins, tres, Nsamples);
 fprintf (candlist, "# Zero-padding the data... %s\n", is_zero_padding ? "yes" : "no");
 fprintf (candlist, "# Search in period range: %lg - %lg ms\n", plow, phigh);
 fprintf (candlist, "#      %s\n", prespres);
 fprintf (candlist, "#      rebinning folded profiles by %ld in the possible range [%ld - %ld]\n", mbins, 1, (p0low * mbins)/2);
 fprintf (candlist, "#      largest extra rebinning of folded profiles by %ld in the possible range [%ld - %ld]\n", pbins, 1, p0low/2);
 fprintf (candlist, "# %s\n", elapsed_time);
 fprintf (candlist, "#\n");
 fprintf (candlist, "# index   period            SNR     phase      Nphases     DC\n");
 fprintf (candlist, "#          (ms)                     (deg)                 (\%)\n");
 fprintf (candlist, "#--------------------------------------------------------------\n");
 for (unsigned long i=0; i<3*pbins*prange; ++i)
  fprintf (candlist, "%5ld    %lf     %8.3lf    %.2lf    %5ld      %.2lf\n", i, period_table[0][i], period_table[1][i], period_table[2][i], (long)period_table[3][i], ((tres * mbins * period_table[4][i]) / period_table[0][i]) * 0.1);
 fclose (candlist);

 double **sift_table;
 sift_table = new double*[5];
 if (!sift_table) { perror("new sift_table"); exit (1); }
 for (int h=0; h<5; ++h) {
  sift_table[h] = new double[3*pbins*prange]; // 3 times larger because we want to keep 3 max snr from every M 
                                              // and then pbins larger because we are doing the same for the range of pbins (from 1 to pbins)
  if (!sift_table[h]) { perror("new sift_table[h]"); exit (1); }
  memset (sift_table[h], 0, 3 * pbins * prange * sizeof(double));
 }
 // sifting period_table
 int sift_size = sifting (period_table, sift_table, 5, 3 * pbins * prange, 2, tres, mbins);

 // output found candidates (only 10 best candidates)
 printf ("# index   period            SNR     phase      Nphases     DC\n");
 printf ("#          (ms)                     (deg)                 (\%)\n");
 printf ("#--------------------------------------------------------------\n");
 for (unsigned long i=0; i<(sift_size > 10 ? 10 : sift_size); ++i) {
  printf ("%5ld    %lf     %8.3lf    %.2lf    %5ld      %.2lf\n", i, sift_table[0][i], sift_table[1][i], sift_table[2][i], (long)sift_table[3][i], ((tres * mbins * sift_table[4][i]) / sift_table[0][i]) * 0.1);
 }

 // forming the name of the sift file
 memset (outname, 0, 255);
 sprintf (outname, "%s/%s%s", output_dir, (strcmp(filestem, "\0") == 0 ? basename(argv[next]) : filestem), ".ffa.sft");

 // writing all sifted period candidates into the file
 FILE *siftlist;
 if ((siftlist = fopen(outname, is_to_append ? "at" : "wt")) == NULL) { perror("fopen"); exit (1); }
 fprintf (siftlist, "# File: %s\n", argv[next]);
 fprintf (siftlist, "# DM = %lg pc/cm^3\n", dm);
 fprintf (siftlist, "# Decimating by %ld   [used tres: %lg mcs] [used Nsamples: %ld]\n", nbins, tres, Nsamples);
 fprintf (siftlist, "# Zero-padding the data... %s\n", is_zero_padding ? "yes" : "no");
 fprintf (siftlist, "# Search in period range: %lg - %lg ms\n", plow, phigh);
 fprintf (siftlist, "#      %s\n", prespres);
 fprintf (siftlist, "#      rebinning folded profiles by %ld in the possible range [%ld - %ld]\n", mbins, 1, (p0low * mbins)/2);
 fprintf (siftlist, "#      largest extra rebinning of folded profiles by %ld in the possible range [%ld - %ld]\n", pbins, 1, p0low/2);
 fprintf (siftlist, "# %s\n", elapsed_time);
 fprintf (siftlist, "#\n");
 fprintf (siftlist, "# index   period            SNR   phase      Nphases     DC\n");
 fprintf (siftlist, "#          (ms)                   (deg)                 (\%)\n");
 fprintf (siftlist, "#--------------------------------------------------------------\n");
 for (unsigned long i=0; i<sift_size; ++i)
  fprintf (siftlist, "%5ld    %lf     %8.3lf    %.2lf    %5ld      %.2lf\n", i, sift_table[0][i], sift_table[1][i], sift_table[2][i], (long)sift_table[3][i], ((tres * mbins * sift_table[4][i]) / sift_table[0][i]) * 0.1);
 fclose (siftlist);

 delete (elapsed_time);
 delete (prespres);
 delete (timepres);
 delete (filestem);
 delete (outname);
 delete (output_dir);
 for (int h=0; h<5; ++h) delete (period_table[h]);
 delete [] period_table;
 for (int h=0; h<5; ++h) delete (sift_table[h]);
 delete [] sift_table;
 
 for (long i=0; i<3; ++i) free (profiles[i]);
 free (profiles);

 for (long i=0; i<3; ++i) free (snr[i]);
 free (snr);

 for (long i=0; i<3; ++i) free (phase[i]);
 free (phase);
 
 for (long i=0; i<psize; ++i) free (y[i]); 
 free (y);
 for (long i=0; i<psize; ++i) free (z[i]); 
 free (z);

 delete (x);
 return 0;
}
