#include "plot.hpp"

/* graphics part */
void plotting (double **array, unsigned long size, char *outname, double plow, double phigh, char *prespres, char *elapsed_time, int argc, char *argv[], int next, double tres, int nbins, int mbins, int pbins, double dm, bool is_zero_padding, double snrmax_p[3][5], float **profiles, SigprocHeader& sigobj) {

 unsigned long p0max = (unsigned long)snrmax_p[0][1];
 if (p0max < (unsigned long)snrmax_p[1][1]) p0max = (unsigned long)snrmax_p[1][1];
 if (p0max < (unsigned long)snrmax_p[2][1]) p0max = (unsigned long)snrmax_p[2][1];

 pgplot_dev = new char[512];
 if (!pgplot_dev) { perror("new pgplot_dev"); exit (1); }
 memset (pgplot_dev, 0, 512);
 sprintf (pgplot_dev, "%s%s%s", outname, ".ps", "/ps");

 char *xlabel = (char *) alloca (255);
 if (!xlabel) { perror("alloca xlabel"); exit (1); }
 char *ylabel = (char *) alloca (255);
 if (!ylabel) { perror("alloca xlabel"); exit (1); }

 float *x, *y, *xprof;
 x = new float[size];
 y = new float[size];
 xprof = new float[p0max];
 if (!x || !y || !xprof) { perror("new x | y | xprof"); exit (1); }
 memset (x, 0, size * sizeof(float));
 memset (y, 0, size * sizeof(float));
 memset (xprof, 0, p0max * sizeof(float));
 for (unsigned long i=0; i<size; ++i) {
  x[i] = (float)array[0][i];
  y[i] = (float)array[1][i];
 }
 float ymax, ymin = y[0];
 for (unsigned long i=0; i<size; ++i) if (y[i] < ymin) ymin = y[i];
 
 ymax = (float)snrmax_p[0][0];

 // multipliers for y-axis
 int ymult[] = { 1000, 500, 250, 200, 150, 100, 50, 25, 10, 5, 4, 2, 1 };
 // multipliers for x-axis
 int xmult[] = { 10000, 5000, 4000, 2000, 1000, 500, 400, 250, 200, 150, 120, 100, 80, 60, 50, 40, 25, 20, 10, 8, 6, 5, 4, 2, 1 };
 double ydataoffset; // data offset in SNR between plot and axis (10% of the y-range)
 double xdataoffset; // offset between y-boundaries and plot points
 ydataoffset = (ymax - ymin) * 0.1;
 xdataoffset = (x[size-1] - x[0]) * 0.02;
 long x_major_tic = gettics (x[0], x[size-1], xmult, sizeof(xmult)/sizeof(int));
 long y_major_tic = gettics (ymin, ymax, ymult, sizeof(ymult)/sizeof(int));

 // plotting SNR-Periods diagram
 cpgopen (pgplot_dev);
 cpgsch (char_size);
 cpgsvp (0.07, 0.49, 0.65, 0.97);
 cpgswin (x[0] - xdataoffset, x[size-1] + xdataoffset, ymin - ydataoffset, ymax + ydataoffset);
 cpgslw (2);
 cpgbox ("BCNTS", x_major_tic, 4, "BCNTS", y_major_tic, 4);
 sprintf (xlabel, "Periods (ms)");
 sprintf (ylabel, "SNR");
 /* x-label */ cpgmtxt ("B", 2.7, 0.5, 0.5, xlabel);
 /* y-label */ cpgmtxt ("L", 2.7, 0.5, 0.5, ylabel);
 cpgslw (8);
 cpgpt (size, x, y, -1);
 cpgpt1 (snrmax_p[0][3], snrmax_p[0][0], 5);
 cpgpt1 (snrmax_p[1][3], snrmax_p[1][0], 5);
 cpgpt1 (snrmax_p[2][3], snrmax_p[2][0], 5);
 cpgslw (1);



 memset (xlabel, 0, 255);
 memset (ylabel, 0, 255);
 sprintf (xlabel, "Time (ms)");
 sprintf (ylabel, "Flux density (\\gs\\fn)");

 float yprofmax = (float)snrmax_p[0][0], yprofmin = profiles[0][0];
 for (unsigned long i=0; i<(unsigned long)snrmax_p[0][1]; ++i) if (profiles[0][i] < yprofmin) yprofmin = profiles[0][i];
 ydataoffset = (yprofmax - yprofmin) * 0.1;
 for (unsigned long i=0; i<p0max; ++i) xprof[i] = i * tres * mbins * snrmax_p[0][4] * 0.001;  // profiles are in ms
 xdataoffset = (xprof[(unsigned long)snrmax_p[0][1]-1] - xprof[0]) * 0.02;

 y_major_tic = gettics (yprofmin, yprofmax, ymult, sizeof(ymult)/sizeof(int));
 x_major_tic = gettics (xprof[0], xprof[(unsigned long)snrmax_p[0][1]-1], xmult, sizeof(xmult)/sizeof(int));

 // Plotting best candidate profile
 cpgsvp (0.57, 0.99, 0.65, 0.97);
 cpgswin (xprof[0] - xdataoffset, xprof[(unsigned long)snrmax_p[0][1]-1] + xdataoffset, yprofmin - ydataoffset, yprofmax + ydataoffset);
 cpgslw (2);
 cpgbox ("BCNTS", x_major_tic, 4, "BCNTS", y_major_tic, 4);
 /* x-label */ cpgmtxt ("B", 2.7, 0.5, 0.5, xlabel);
 /* y-label */ cpgmtxt ("L", 2.7, 0.5, 0.5, ylabel);
 cpgslw (2);
 cpgline ((unsigned long)snrmax_p[0][1], xprof, profiles[0]);

 yprofmax = (float)snrmax_p[1][0]; 
 yprofmin = profiles[1][0];
 for (unsigned long i=0; i<(unsigned long)snrmax_p[1][1]; ++i) if (profiles[1][i] < yprofmin) yprofmin = profiles[1][i];
 ydataoffset = (yprofmax - yprofmin) * 0.1;
 for (unsigned long i=0; i<p0max; ++i) xprof[i] = i * tres * mbins * snrmax_p[1][4] * 0.001;  // profiles are in ms
 xdataoffset = (xprof[(unsigned long)snrmax_p[1][1]-1] - xprof[0]) * 0.02;
 y_major_tic = gettics (yprofmin, yprofmax, ymult, sizeof(ymult)/sizeof(int));
 x_major_tic = gettics (xprof[0], xprof[(unsigned long)snrmax_p[1][1]-1], xmult, sizeof(xmult)/sizeof(int));

 // Plotting second candidate profile
 cpgsvp (0.07, 0.49, 0.25, 0.57);
 cpgswin (xprof[0] - xdataoffset, xprof[(unsigned long)snrmax_p[1][1]-1] + xdataoffset, yprofmin - ydataoffset, yprofmax + ydataoffset);
 cpgslw (2);
 cpgbox ("BCNTS", x_major_tic, 4, "BCNTS", y_major_tic, 4);
 /* x-label */ cpgmtxt ("B", 2.7, 0.5, 0.5, xlabel);
 /* y-label */ cpgmtxt ("L", 2.7, 0.5, 0.5, ylabel);
 cpgslw (2);
 cpgline ((unsigned long)snrmax_p[1][1], xprof, profiles[1]);

 yprofmax = (float)snrmax_p[2][0]; 
 yprofmin = profiles[2][0];
 for (unsigned long i=0; i<(unsigned long)snrmax_p[2][1]; ++i) if (profiles[2][i] < yprofmin) yprofmin = profiles[2][i];
 ydataoffset = (yprofmax - yprofmin) * 0.1;
 for (unsigned long i=0; i<p0max; ++i) xprof[i] = i * tres * mbins * snrmax_p[2][4] * 0.001;  // profiles are in ms
 xdataoffset = (xprof[(unsigned long)snrmax_p[2][1]-1] - xprof[0]) * 0.02;
 y_major_tic = gettics (yprofmin, yprofmax, ymult, sizeof(ymult)/sizeof(int));
 x_major_tic = gettics (xprof[0], xprof[(unsigned long)snrmax_p[2][1]-1], xmult, sizeof(xmult)/sizeof(int));

 // Plotting third candidate profile
 cpgsvp (0.57, 0.99, 0.25, 0.57);
 cpgswin (xprof[0] - xdataoffset, xprof[(unsigned long)snrmax_p[2][1]-1] + xdataoffset, yprofmin - ydataoffset, yprofmax + ydataoffset);
 cpgslw (2);
 cpgbox ("BCNTS", x_major_tic, 4, "BCNTS", y_major_tic, 4);
 /* x-label */ cpgmtxt ("B", 2.7, 0.5, 0.5, xlabel);
 /* y-label */ cpgmtxt ("L", 2.7, 0.5, 0.5, ylabel);
 cpgslw (2);
 cpgline ((unsigned long)snrmax_p[2][1], xprof, profiles[2]);

 /* printing legend */
 cpgsvp (0.07, 0.99, 0.25, 0.97);
 cpgsch (0.8*char_size);

 char *cur;
 cur = new char[255];
 if (!cur) { perror("new cur"); exit (1); }

 memset (cur, 0, 255);
 if (sigobj.is_initialized()) {
  char *telescope, *source, *date;
  telescope = new char[80];
  source = new char[80];
  date = new char[11];
  if (!telescope || !source || !date) { perror ("new telescope | source | date"); exit (1); }
  memset (telescope, 0, 80);
  memset (source, 0, 80);
  memset (date, 0, 11);
  telescope = sigobj.get_telescope (telescope);
  source = sigobj.get_source (source);
  date = sigobj.get_date (date);

  sprintf (cur, "%s   Source: %s   MJD: %lf  Date: %s", telescope, source, sigobj.get_mjd(), date); 

  delete (telescope);
  delete (source);
  delete (date);
 } else sprintf (cur, "Telescope: %s  Source: %s  MJD: %s  Date: %s", "", "", "", "");
 cpgmtxt ("B", 5., 0., 0., cur);

 memset (cur, 0, 255);
 if (sigobj.is_initialized()) sprintf (cur, "Data file: %s  Observation length: %.1f s", argv[next], sigobj.get_tobs());
  else sprintf (cur, "Data file: %s  Observation length: %s", argv[next], "");
 cpgmtxt ("B", 6., 0., 0., cur);

 memset (cur, 0, 255);
 sprintf (cur, "Original tres: %lg \\gm\\fns   Decim by %d   Used tres: %lg \\gm\\fns   DM: %lg pc/cm\\u3\\d   Zero-pad: %s   mbins: %d   pbins: %d", tres/nbins, nbins, tres, dm, is_zero_padding ? "yes" : "no", mbins, pbins);
 cpgmtxt ("B", 7., 0., 0., cur);

 memset (cur, 0, 255);
 sprintf (cur, "FFA Search:  period range: %lg - %lg ms  %s", plow, phigh, prespres);
 cpgmtxt ("B", 8., 0., 0., cur);

 memset (cur, 0, 255);
 sprintf (cur, "Best period candidates: %lf (snr = %.3lf)   %lf (%.3lf)   %lf (%.3lf)", snrmax_p[0][3], snrmax_p[0][0], snrmax_p[1][3], snrmax_p[1][0], snrmax_p[2][3], snrmax_p[2][0]);
 cpgmtxt ("B", 9., 0., 0., cur);

 cpgsch (0.7*char_size);
 memset (cur, 0, 255);
 sprintf (cur, "Cmd options: ");
 for (int i=0; i<next; i++) sprintf (cur, "%s %s", cur, argv[i]);
 cpgmtxt ("B", 12.0, 0., 0., cur);

 // print name of the output device (file)
 cpgsch (0.6*char_size);
 cpgmtxt ("B", 16., 0., 0., elapsed_time);
 memset (cur, 0, 255);
 sprintf (cur, "%s%s", basename (outname), ".ps");
 cpgmtxt ("B", 16., 1.01, 1., cur);
 cpgclos();

 delete (pgplot_dev);
 delete (x);
 delete (y);
 delete (cur);
}


// get tics major interval
long gettics (float min, float max, int mult[], int sz) {
 long beg = (long)floor((max - min)/5.);
 long end = (long)floor((max - min)/4.);
 for (int j=0; j<sz; ++j)
  for (long i=beg; i<=end; ++i) {
   if (((long)((float)i/(float)mult[j]))*mult[j] == i) return i;
  }
 return end;
}
