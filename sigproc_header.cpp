/* 
 * I took almost all the code here from the file read_header.c from
 * sigproc package written by Duncan Lorimer
 * to get wrapped into the SigprocHeader class
 * some functions are also taken from aliases.c and slalib.c
 * and nsamples.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>

#include "sigproc_header.hpp"
#include "slamac.h"

/* I have moved this from strings_equal.c file in sigproc package. Vlad */
int SigprocHeader::strings_equal (char *string1, char *string2) { /* includefile */
 if (!strcmp(string1, string2)) return 1;
   else return 0;
}

/* read a string from the input which looks like nchars-char[1-nchars] */
void SigprocHeader::get_string(int f, int *nbytes, char string[]) { /* includefile */
 int nchar;
 strcpy (string, "ERROR");
 read (f, &nchar, sizeof(int));
 struct stat info;
 if (fstat(f, &info)) { perror("stat"); exit (1); }
 if (nchar >= info.st_size) return;
 if (nchar > 80 || nchar < 1) return;
 *nbytes=sizeof(int);
 read(f, string, nchar);
 string[nchar]='\0';
 *nbytes+=nchar;
}

/* constructor */
SigprocHeader::SigprocHeader (char *inputfile) {
 is_sigproc = false;
 initialize (inputfile);
}

/* this function initialize the fields of the header from the input file */
/* attempt to read in the general header info from a pulsar data file */
void SigprocHeader::initialize (char *inputfile) {

  char string[80], message[80];
  int itmp, nbytes, totalbytes, expecting_rawdatafile = 0, expecting_source_name = 0; 
  int expecting_frequency_table = 0, channel_index;

  int f; // file descriptor
  if ((f = open (inputfile, O_RDONLY | O_LARGEFILE)) == -1) { perror("open"); exit(1); };
  lseek (f, 0, SEEK_SET);

  struct stat inf;
  if (fstat(f, &inf)) { perror("stat"); exit (1); }
  off_t filesize = inf.st_size;

  /* try to read in the first line of the header */
  get_string (f, &nbytes, string);
  if (!strings_equal(string, "HEADER_START")) {
	/* the data file is not in standard format, rewind and return */
	lseek (f, 0, SEEK_SET);
	return;
  }
  /* store total number of bytes read so far */
  totalbytes=nbytes;

  /* loop over and read remaining header lines until HEADER_END reached */
  while (1) {
    get_string (f, &nbytes, string);
    if (strings_equal (string, "HEADER_END")) break;
    totalbytes+=nbytes;
    if (strings_equal (string, "rawdatafile")) {
      expecting_rawdatafile=1;
    } else if (strings_equal (string, "source_name")) {
      expecting_source_name=1;
    } else if (strings_equal (string, "FREQUENCY_START")) {
      expecting_frequency_table=1;
      channel_index=0;
    } else if (strings_equal (string, "FREQUENCY_END")) {
      expecting_frequency_table=0;
    } else if (strings_equal (string, "az_start")) {
      read(f, &az_start, sizeof(az_start));
      totalbytes+=sizeof(az_start);
    } else if (strings_equal (string, "za_start")) {
      read(f, &za_start, sizeof(za_start));
      totalbytes+=sizeof(za_start);
    } else if (strings_equal (string, "src_raj")) {
      read(f, &src_raj, sizeof(src_raj));
      totalbytes+=sizeof(src_raj);
    } else if (strings_equal (string, "src_dej")) {
      read(f, &src_dej, sizeof(src_dej));
      totalbytes+=sizeof(src_dej);
    } else if (strings_equal (string, "tstart")) {
      read(f, &tstart, sizeof(tstart));
      totalbytes+=sizeof(tstart);
    } else if (strings_equal (string, "tsamp")) {
      read(f, &tsamp, sizeof(tsamp));
      totalbytes+=sizeof(tsamp);
    } else if (strings_equal (string, "period")) {
      read(f, &period, sizeof(period));
      totalbytes+=sizeof(period);
    } else if (strings_equal (string, "fch1")) {
      read(f, &fch1, sizeof(fch1));
      totalbytes+=sizeof(fch1);
    } else if (strings_equal (string, "fchannel")) {
      read(f, &frequency_table[channel_index++], sizeof(double));
      totalbytes+=sizeof(double);
      fch1=foff=0.0; /* set to 0.0 to signify that a table is in use */
    } else if (strings_equal (string, "foff")) {
      read(f, &foff, sizeof(foff));
      totalbytes+=sizeof(foff);
    } else if (strings_equal (string, "nchans")) {
      read(f, &nchans, sizeof(nchans));
      totalbytes+=sizeof(nchans);
    } else if (strings_equal (string, "telescope_id")) {
      read(f, &telescope_id, sizeof(telescope_id));
      totalbytes+=sizeof(telescope_id);
    } else if (strings_equal (string, "machine_id")) {
      read(f, &machine_id, sizeof(machine_id));
      totalbytes+=sizeof(machine_id);
    } else if (strings_equal (string, "data_type")) {
      read(f, &data_type, sizeof(data_type));
      totalbytes+=sizeof(data_type);
    } else if (strings_equal (string, "ibeam")) {
      read(f, &ibeam, sizeof(ibeam));
      totalbytes+=sizeof(ibeam);
    } else if (strings_equal (string, "nbeams")) {
      read(f, &nbeams, sizeof(nbeams));
      totalbytes+=sizeof(nbeams);
    } else if (strings_equal (string, "nbits")) {
      read(f, &nbits, sizeof(nbits));
      totalbytes+=sizeof(nbits);
    } else if (strings_equal (string, "barycentric")) {
      read(f, &barycentric, sizeof(barycentric));
      totalbytes+=sizeof(barycentric);
    } else if (strings_equal (string, "pulsarcentric")) {
      read(f, &pulsarcentric, sizeof(pulsarcentric));
      totalbytes+=sizeof(pulsarcentric);
    } else if (strings_equal (string, "nbins")) {
      read(f, &nbins, sizeof(nbins));
      totalbytes+=sizeof(nbins);
    } else if (strings_equal (string, "nsamples")) {
      /* read this one only for backwards compatibility */
      read(f, &itmp, sizeof(itmp));
      totalbytes+=sizeof(itmp);
    } else if (strings_equal (string, "nifs")) {
      read(f, &nifs, sizeof(nifs));
      totalbytes+=sizeof(nifs);
    } else if (strings_equal (string, "npuls")) {
      read(f, &npuls, sizeof(npuls));
      totalbytes+=sizeof(npuls);
    } else if (strings_equal (string, "refdm")) {
      read(f, &refdm, sizeof(refdm));
      totalbytes+=sizeof(refdm);
    } else if (expecting_rawdatafile) {
      strcpy(rawdatafile, string);
      expecting_rawdatafile=0;
    } else if (expecting_source_name) {
      strcpy(source_name, string);
      expecting_source_name=0;
    } else {
      sprintf(message,"read_header - unknown parameter: %s\n", string);
      fprintf(stderr, "ERROR: %s\n", message);
      return;
    } 
  } 
  // close the file
  close (f);

  /* add on last header string */
  totalbytes+=nbytes;

  /* set total number of bytes read */
  header_size = totalbytes;

  // calc the number of samples
  nsamp = nsamples (filesize);
  tobs = nsamp * tsamp;
  calc_date (tstart); // calculate the date from tstart, setting date[10] array

  // reading of the header is successful
  is_sigproc = true;
}

// calculate the number of samples in the file
// code was taken from nsamples.c from sigproc package
off_t SigprocHeader::nsamples (off_t filesize) {
  off_t datasize, numsamps;
  datasize=filesize - header_size;
  numsamps= (off_t) ((long double) (datasize) / (((long double) nbits) / 8.0) / (long double) nifs / (long double) nchans);
  return numsamps;
}

// get telescope name
// code taken from aliases.c from sigproc package
char* SigprocHeader::get_telescope (char *dest) {
 switch (telescope_id) {
   case 0: strcpy (dest, "Fake");
	   break;
   case 1: strcpy (dest, "Arecibo");
           break;
   case 2: strcpy (dest, "Ooty");
           break;
   case 3: strcpy (dest, "Nancay");
           break;
   case 4: strcpy (dest, "Parkes");
           break;
   case 5: strcpy (dest, "Jodrell");
           break;
   case 6: strcpy (dest, "GBT");
           break;
   case 7: strcpy (dest, "GMRT");
           break;
   case 8: strcpy (dest, "Effelsberg");
           break;
   default: strcpy (dest, "???????");
            break;
 } // switch
 return dest;
}

// get date
char* SigprocHeader::get_date (char *dest) {
 strcpy (dest, date);
 return dest;
}

// get source name
char* SigprocHeader::get_source (char *dest) {
 strcpy (dest, source_name);
 return dest;
}

// to set the date from given MJD
void SigprocHeader::calc_date (double mjd) {
  int iymdf[4], j;
  slaDjcal (1, mjd, iymdf, &j);
  if (j == 0) sprintf (date, "%04d/%02d/%02d", iymdf[0], iymdf[1], iymdf[2]);
   else sprintf (date, "%s/%s/%s", "????", "??", "??");
}

// to calculate the date from the given MJD
// this code is taken from slalib.c from sigproc package
void SigprocHeader::slaDjcal (int ndp, double djm, int iymdf[4], int *j) {
/*
**  - - - - - - - - -
**   s l a D j c a l
**  - - - - - - - - -
**
**  Modified Julian Date to Gregorian calendar, expressed
**  in a form convenient for formatting messages (namely
**  rounded to a specified precision, and with the fields
**  stored in a single array).
**
**  Given:
**     ndp      int       number of decimal places of days in fraction
**     djm      double    Modified Julian Date (JD-2400000.5)
**
**  Returned:
**     iymdf    int[4]    year, month, day, fraction in Gregorian calendar
**     *j       int       status:  nonzero = out of range
**
**  Any date after 4701BC March 1 is accepted.
**
**  Large ndp values risk internal overflows.  It is typically safe
**  to use up to ndp=4.
**
**  The algorithm is derived from that of Hatcher 1984 (QJRAS 25, 53-55).
**
**  Defined in slamac.h:  dmod
**
**  Last revision:   17 August 1999
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
   double fd, df, f, d;
   long jd, n4, nd10;

/* Validate */
   if ( ( djm <= -2395520.0 ) || ( djm >= 1.0e9 ) ) {
      *j = - 1;
      return;
   } else {

   /* Denominator of fraction */
      fd = pow ( 10.0, (double) gmax ( ndp, 0 ) );
      fd = dnint ( fd );

   /* Round date and express in units of fraction */
      df = djm * fd;
      df = dnint ( df );

   /* Separate day and fraction */
      f = dmod ( df, fd );
      if ( f < 0.0 ) f += fd;
      d = ( df - f ) / fd;

   /* Express day in Gregorian calendar */
      jd = (long) dnint ( d ) + 2400001L;
      n4 = 4L * ( jd + ( ( 2L * ( ( 4L * jd - 17918L ) / 146097L)
                                       * 3L ) / 4L + 1L ) / 2L - 37L );
      nd10 = 10L * ( ( ( n4 - 237L ) % 1461L ) / 4L ) + 5L;
      iymdf[0] = (int) ( ( n4 / 1461L ) - 4712L );
      iymdf[1] = (int) ( ( ( nd10 / 306L + 2L ) % 12L ) + 1L );
      iymdf[2] = (int) ( ( nd10 % 306L ) / 10L + 1L );
      iymdf[3] = (int) dnint ( f );
      *j = 0;
   }
}
