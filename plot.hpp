#ifndef _PLOT_HPP
#define _PLOT_HPP

#include <stdlib.h>

#include "cpgplot.h"
#include "sigproc_header.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <alloca.h>
#include <libgen.h>

extern char *pgplot_dev; // device for pgplot
extern float char_size;  // character size

/* graphics part */
void plotting (double **array, unsigned long size, char *outname, double plow, double phigh, char *prespres, char *elapsed_time, int argc, char *argv[], int next, double tres, int nbins, int mbins, int pbins, double dm, bool is_zero_padding, double snrmax_p[3][5], float **profiles, SigprocHeader& sigobj);

// get tics major interval
long gettics (float min, float max, int mult[], int sz);

#endif //#ifndef _PLOT_HPP
