#include "cmd.hpp"

/* help */ 
void help (char *s) {
 printf("\nPROGRAMM: %s\n"
	"Run a FFA search in time series data\n\n"
        "[USAGE]: %s [options] <.tim|.dat>\n\n"
        "[OPTIONS]:\n"
	"--plow <period>     - the lowest period to search for in ms (default = 1000)\n"
	"--phigh <period>    - the highest period to search for in ms (default = 10000)\n"
	"-t, --tres <value>  - time resolution in mcs (default = 81.92 or read from the sigproc header)\n"
	"--dm <value>        - dispersion measure in pc/cm^3 (default = 0 or read from the sigproc header)\n"
	"-z, --zero          - to zero-pad the data to the nearest power of 2\n"
	"--header <bytes>    - size of the header in bytes (default = 0)\n"
	"-l, --left <pos>    - _pos_ number of samples to skip in the input file, default = 0\n"
	"-w, --window <size> - _size_ number of samples to read from the file, default - all samples\n"
	"--periodres <value> - period resolution _value_ in mcs, default - it is auto calculated as tres/M, where\n"
	"                      M is the nearest down of power of 2 of the ratio Nsamples/period (in counts).\n"
	"                      Since FFA requires M not to be less than 2, period resolution should not be larger than half of tres!\n"
	"                      This option fixes the only resolution in periods in the whole range of periods,\n"
	"                      Nsamples will be auto calculated for every current processed period\n"
	"                      This option is not implemented yet completely, only first block will be read\n"
	"--mres <value>      - _value_ of M, which should be power of 2, see --periodres\n"
	"                      If used with --periodres, the lower value of M will be used\n"
	"                      This option is not implemented yet completely, only first block will be read\n"
	"-b, --block <size>  - maximum _size_ number of bytes to be read at once from the input file, default - 400000000\n"
	"-s, --sigproc       - consider the tim-file written in sigroc format\n"
	"-n, --nbins <value> - to downsample the raw-data time series by factor _value_ (default = 1)\n"
	"-m, --mbins <value> - to rebin the folded profiles by factor _value_ in FFA (default = 1)\n"
	"                      this option determines the lowest duty cycle of the pulsar\n"
	"-p, --pbins <value> - largest extra rebinning of the folded profiles by factor _value_\n"
	"                      after the FFA before looking for candidates (default = 1)\n"
	"                      this option determines the largest duty cycle to check. Program will search\n"
	"                      in the range of rebin factors from 1 to _value_ with the step 1\n"
	"-d, --dir <dir>     - output directory for candidate list file, default - current directory\n"
	"--filestem <string> - filestem used to form the name of output files, default - full name of input file\n"
	"--append            - to append output candidate list file, default - to overwrite\n"
	"--verbose           - to show percentage of processed trial periods\n"
	"-h, --help          - help\n\n", basename(s), basename(s));
 exit(0);
}

/* parse command line */
int parse_command_line (int argc, char *argv[]) {
  int op, tm;
  struct option long_options[] = { {"help", no_argument, 0, 'h'},
  	                           {"plow", required_argument, 0, 10},
  	                           {"phigh", required_argument, 0, 11},
  	                           {"tres", required_argument, 0, 't'},
  	                           {"dm", required_argument, 0, 12},
				   {"zero", no_argument, 0, 'z'},
  	                           {"header", required_argument, 0, 13},
				   {"sigproc", no_argument, 0, 's'},
  	                           {"nbins", required_argument, 0, 'n'},
  	                           {"mbins", required_argument, 0, 'm'},
  	                           {"pbins", required_argument, 0, 'p'},
				   {"verbose", no_argument, 0, 14},
  	                           {"dir", required_argument, 0, 'd'},
  	                           {"filestem", required_argument, 0, 15},
				   {"append", no_argument, 0, 16},
  	                           {"left", required_argument, 0, 'l'},
  	                           {"window", required_argument, 0, 'w'},
  	                           {"block", required_argument, 0, 'b'},
  	                           {"periodres", required_argument, 0, 17},
  	                           {"mres", required_argument, 0, 18},
				   {0, 0, 0, 0}
                                  };

  optind = 0; optarg = 0;
  while((op = getopt_long(argc, argv, "t:zsn:m:d:l:w:b:p:h", long_options, 0)) != EOF)
    switch(op){
      
      case 'h':
        help (argv[0]);
      break;

      case 10:
     	plow = atof (optarg); 
      break;
      
      case 11:
     	phigh = atof (optarg); 
      break;
      
      case 't':
     	tres = atof (optarg); 
      break;
      
      case 12:
     	dm = atof (optarg); 
      break;
      
      case 'z':
     	is_zero_padding = true;
      break;
      
      case 13:
     	if (atoi (optarg) < 0) header = 0; 
	 else header = atoi (optarg);
      break;
      
      case 's':
     	is_sigproc = true;
      break;
      
      case 'n':
     	nbins = atoi (optarg); 
	if (nbins < 1) nbins = 1;
      break;

      case 'm':
     	mbins = atoi (optarg); 
	if (mbins < 1) mbins = 1;
      break;

      case 'p':
     	pbins = atoi (optarg); 
	if (pbins < 1) pbins = 1;
      break;

      case 14:
     	is_verbose = true;
      break;

      case 'd':
        if (output_dir) { output_dir = 0; free (output_dir); }
        if (!(output_dir = (char *)calloc((int)strlen(optarg), sizeof(char *)))) { perror("new output_dir"); exit (1); }
        memcpy (output_dir, optarg, (int)strlen(optarg));
      break;

      case 15:
        if (filestem) { filestem = 0; free (filestem); }
        if (!(filestem = (char *)calloc((int)strlen(optarg), sizeof(char *)))) { perror("new filestem"); exit (1); }
        memcpy (filestem, optarg, (int)strlen(optarg));
      break;

      case 16:
     	is_to_append = true;
      break;

      case 'l':
     	left_edge = atol (optarg); 
	if (left_edge < 0) left_edge = 0;
      break;

      case 'w':
      	Nsamples = atol (optarg);
	if (Nsamples < 1) Nsamples = 0;
      break;
      
      case 'b':
	if (atol(optarg) >= 1) block = atol (optarg);
      break;

      case 17:
      	period_res = atof (optarg);
	is_period_res = true;
      break;

      case 18:
      	Mcmd = atol (optarg);
	is_Mcmd = true;
      break;

      case '?':
        help (argv[0]);
      break;
    }
  return(optind);
}
