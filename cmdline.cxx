#include <stdio.h>
#include <unistd.h>
#include <libgen.h>
#include <getopt.h>
#include "dragnet.h"

// Usage help
void usage(char *prg) {
 printf("USAGE: %s [OPTIONS] <file>\n\n", basename(prg));
 printf("OPTIONS:\n");
 printf(" -h, --help                    This help\n");
 printf(" -q, --quiet                   Quiet; no information to screen\n");
 printf(" -f, --format <string>         Input data format. Possible values are: sigproc, hdf5. Default - sigproc\n");
 printf(" -b, --blocksize <int>         Maximum number of samples to read at once. Default - whole file\n");
 printf(" -o, --output <prefix>         Output prefix. Default - test\n");
 printf(" -D, --gpu <id>                GPU device number to use. Default - 0\n");
 printf(" -r, --range <startDM,endDM>   Start and end values of DM range to dedisperse. Default - 0,50\n");
 printf(" -s, --dmstep <DMstep>         Linear DM step to use. If not provided, optimal DM trials are computed\n");
 printf(" -w, --width <pulseWidth>      Expected intrinsic pulse width for optimal DM trials. Default - 4.0us\n");
 printf(" -t, --tolerance <tolerance>   Smearing tolerance factor between DM trials. Default - 1.25\n");
 printf(" -m, --mask <maskfile>         Maskfile *_rfifind.mask to apply\n");
 printf(" -c, --clip <threshold>        Clip all samples above <threshold> in sigmas. Only can be used with --mask. Default - 0 (no clipping)\n");
 printf(" -z, --zapchan <list>          Comma separated string (no spaces!) of channels to explicitly remove from analysis.\n");
 printf("                               Ranges are specified by min:max[:step]. The lower channel number, the lower the frequency\n");
 exit(0);
}

// parsing command line arguments
// return value is the index of the first non-option argument (first input file)
int parse_cmdline(int argc, char *argv[], cmdline* cmd) {
  int op;
  struct option long_options[] = { {"help", no_argument, 0, 'h'},
                                   {"quiet", no_argument, 0, 'q'},
                                   {"format", required_argument, 0, 'f'},
                                   {"output", required_argument, 0, 'o'},
                                   {"gpu", required_argument, 0, 'D'},
                                   {"range", required_argument, 0, 'r'},
                                   {"dmstep", required_argument, 0, 's'},
                                   {"width", required_argument, 0, 'w'},
                                   {"tolerance", required_argument, 0, 't'},
                                   {"blocksize", required_argument, 0, 'b'},
                                   {"mask", required_argument, 0, 'm'},
                                   {"clip", required_argument, 0, 'c'},
                                   {"zapchan", required_argument, 0, 'z'},
                                   {0, 0, 0, 0}
                                  };

  while((op = getopt_long(argc, argv, "hqf:o:D:r:s:w:t:b:m:c:z:", long_options, 0)) != EOF)
    switch(op){
      case 'h':
        usage(argv[0]);
      break;

      case 'q':
	cmd->verbose = 0;
      break;

      case 'f':
        strcpy(cmd->format, optarg);
      break;

      case 'o':
        strcpy(cmd->prefix, optarg);
      break;

      case 'D':
        cmd->device_id = atoi(optarg);
      break;

      case 'r':
        sscanf(optarg,"%f,%f",&(cmd->dm_start), &(cmd->dm_end));
      break;

      case 's':
        cmd->dm_step = atof(optarg);
      break;

      case 'w':
        cmd->pulse_width = atof(optarg);
      break;

      case 't':
        cmd->dm_tol = atof(optarg);
      break;

      case 'b':
        cmd->blocksize = strtoull(optarg, NULL, 10);
      break;

      case 'm':
        strcpy(cmd->maskfile, optarg);
      break;

      case 'c':
        cmd->clip_sigma = atof(optarg);
      break;

      case 'z':
        strcpy(cmd->zapchan, optarg);
      break;

      case '?':
        usage(argv[0]);
      break;
    }
  return(optind);
 }
