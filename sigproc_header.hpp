#ifndef _SIGPROC_HEADER_HPP
#define _SIGPROC_HEADER_HPP

#include <stdio.h>

/*
 * This class describes the fields of the header of sigproc package
 * written by Duncan Lorimer
 * I have taken almost all the code from header.h and read_header.c files
 * of this package. Vlad
 */
class SigprocHeader {

int header_size; // size of the header in bytes
bool is_sigproc; // the flag of the success of reading the header, true - success
char date[10];   // date of the observations in format YYYY/MM/DD
float tobs;   // length of observations in seconds
off_t nsamp;  // number of samples

/* global variables describing the data */
char rawdatafile[80]; 
char source_name[80];
int machine_id; 
int telescope_id; 
int data_type; 
int nbeams; 
int ibeam;
int nbits; 

/* these two added Aug 20, 2004 DRL */
int barycentric;
int pulsarcentric;

int nchans; 
int nbins;
int nifs; 
long int npuls; /* added for binary pulse profile format */

double az_start;
double za_start;
double src_raj;
double src_dej;
double tstart;
double tsamp;
double period;
double fch1;
double foff;
double refdm;

/* added frequency table for use with non-contiguous data */
double frequency_table[4096]; /* note limited number of channels */

  int strings_equal (char *string1, char *string2); // just to compare strings
  void get_string(int f, int *nbytes, char string[]); // read string from the file
  // calculate the number of samples in the file
  // // code was taken from nsamples.c from sigproc package
  off_t nsamples (off_t filesize);
  // to calculate the date from the given MJD
  // this code is taken from slalib.c from sigproc package
  void slaDjcal (int ndp, double djm, int iymdf[4], int *j);
  // to set the date from given MJD
  void calc_date (double mjd);

  public:
   SigprocHeader () { is_sigproc = false; }; // empty constructor
   SigprocHeader (char *inputfile);
   virtual ~SigprocHeader () {}

   /* this function initialize the fields of the header from the input file */
   void initialize (char *inputfile);

   // get the size of the header
   int get_header_size () { return header_size; }

   // check if sigproc header was read successfully
   bool is_initialized () { return is_sigproc; }
   // get reference DM
   double get_refdm () { return refdm; }
   // get sampling time
   double get_tsamp () { return tsamp; }
   // get date
   char* get_date (char *dest);
   // get length of observations
   float get_tobs () { return tobs; }
   // get source name
   char* get_source (char *dest);
   // get telescope name
   char* get_telescope (char *dest);
   // get mjd of the first sample
   double get_mjd () { return tstart; }
};

#endif //#ifndef _SIGPROC_HEADER_HPP
