#ifndef _FFASEARCH_HPP
#define _FFASEARCH_HPP

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <libgen.h>

#if TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  include <time.h>
# endif
#endif

#include "cmd.hpp"
#include "dosum.hpp"
#include "dosnr.hpp"
#include "rebin.hpp"
#include "plot.hpp"
#include "sigproc_header.hpp"

// define maximum number of samples to read from file at once
// in bytes, default - 100 mln. bytes
// it's of about 70% of total memory in 2GB systems
#define MAXBLOCKSIZE 400000000

double tres = 81.92;          // time resolution in mcs (to be read from the header)
unsigned int header = 0;      // size of the header (in bytes)
double plow = 1000.;          // lowest period in the range to search (in ms)
double phigh = 10000.;        // highest period in the range to be search (in ms)
int nbins = 1;                // number of samples to add (rebin factor for raw-data)
int mbins = 1;                // number of samples to add in folded profiles (for searching)
int pbins = 1;                // number of samples to add in folded profiles after the FFA before looking for candidates
double dm = 0.;               // dispersion measure in pc/cm^3
bool is_sigproc = false;      // if true, tim-file is written in sigproc format
bool is_zero_padding = false; // if true, use zero-padding to the nearest power of 2
bool is_verbose = false;      // if true, number of processed (in %) trial periods will be shown
bool is_to_append = false;    // if true, output candidate list file will be appended to exist one
                              // otherwise, it will be truncated
double wpratio = 0.03;        // the ratio of the pulse width to the period. I assume that in general it is of about 3%
                              // It uses to estimate the effective pulse width and to suggest the rebinning factor
			      // And also it is used to sift the candidates, and if their phases fall into the interval of 3% 
			      // of the period than these candidates are considered to be the same

off_t filesamples;   // number of samples in tim-file (to be read from the header)
off_t left_edge = 0; // left pos in samples to read from the file, default = 0
off_t Nsamples = 0;  // number of samples to read from the file, default - all file

off_t block = MAXBLOCKSIZE; // size of the block (in bytes) to be read at once from the file

char *output_dir;  // output dir
char *filestem;    // filestem

double period_res; // resolution in period (it depends from the base period and number of samples)
bool is_period_res = false;  // it true, period_res was set in command line and it will be used

unsigned long Mcmd; // number of M in command line, if used together with --periodres than lower value of M will be used
bool is_Mcmd = false;  // if M is set in command line

float **y; // 2-D array for FFA search
float **z; // backup 2-D array

SigprocHeader sigobj; // object for sigproc header

/*--- plot.cpp ---*/
char *pgplot_dev;     // device for pgplot
float char_size = 1.; // character size

#endif //#ifndef _FFASEARCH_HPP
