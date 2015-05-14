#ifndef _CMD_HPP
#define _CMD_HPP

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

extern double tres;          // time resolution in mcs (to be read from the header)
extern unsigned int header;  // size of the header (in bytes)
extern double plow;          // lowest period in the range to search (in ms)
extern double phigh;         // highest period in the range to be search (in ms)
extern int nbins;            // number of samples to add (rebin factor for raw-data)
extern int mbins;            // number of samples to add in folded profiles (for searching)
extern int pbins;            // number of samples to add in folded profiles after FFA before looking for candidates
extern double dm;            // dispersion measure
extern bool is_sigproc;      // if true, tim-file is written in sigproc format
extern bool is_zero_padding; // if true, use zero-padding to the nearest power of 2
extern bool is_verbose;      // if true, number of processed (in %) trial periods will be shown 
extern char* output_dir;     // output directory
extern char* filestem;       // filestem
extern bool is_to_append;    // if true, output candidate list file will be appended to exist one
                             // otherwise, it will be truncated
extern off_t left_edge;      // left pos in samples to read from the file, default = 0
extern off_t Nsamples;       // number of samples to read from the file, default - all file

extern off_t block;          // size of the block (in bytes) to be read at once from the file

extern double period_res;    // resolution in period (it depends from the base period and number of samples)
extern bool is_period_res;   // it true, period_res was set in command line and it will be used

extern unsigned long Mcmd;   // number of M in command line, if used together with --periodres than lower value of M will be used
extern bool is_Mcmd;         // if M is set in command line

/* help */ 
void help (char *s);

/* parse command line */
int parse_command_line (int argc, char *argv[]);

#endif //#ifndef _CMD_HPP
