#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <dedisp.h>

#ifndef _HEADER
#define _HEADER

typedef struct {
  int64_t headersize,buffersize;
  unsigned int nchan,nsamp,nbit,nif;
  int machine_id,telescope_id,nbeam,ibeam,sumif;
  double tstart,tsamp,fch1,foff;
  double src_raj,src_dej,az_start,za_start;
  char source_name[80],ifstream[8],inpfile[80];
} header;

#endif

#ifndef _CMDLINE
#define _CMDLINE

typedef struct {
 int device_id, verbose;
 unsigned long long blocksize;
 char prefix[128], format[16];
 dedisp_float dm_start, dm_end, dm_step, pulse_width, dm_tol;
} cmdline;

#endif

/*-- inf.cxx --*/
char *backend_name(int machine_id);
char *telescope_name(int telescope_id);
void writeinf(header* h, char *outstem, float dm, int64_t shift_back);

/*-- input.cxx --*/
int raw_open(char *filename, char *format, header *h, int verbose, void*& fd); // open input file
int64_t raw_read(int64_t nsamples, int64_t shift_back, header* h, void*& out, void*& fd); // read data
void raw_close(void*& fd); // close input file

// Usage help
void usage(char *prg);
// parsing command line arguments
int parse_cmdline(int argc, char *argv[], cmdline* cmd);
