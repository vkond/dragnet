#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <dedisp.h>

typedef struct {
  int64_t headersize,buffersize;
  int nchan,nsamp,nbit,nif;
  int machine_id,telescope_id,nbeam,ibeam,sumif;
  double tstart,tsamp,fch1,foff;
  double src_raj,src_dej,az_start,za_start;
  char source_name[80],ifstream[8],inpfile[80];
} header;

typedef struct {
 int device_id, verbose;
 char prefix[128], format[16];
 dedisp_float dm_start, dm_end, dm_step, pulse_width, dm_tol;
} cmdline;

// reading Sigproc header
int read_sigproc_header(FILE* file, header* h);
static char *backend_name(int machine_id);
static char *telescope_name(int telescope_id);
void writeinf(header* h, char *outstem, float dm);

// opens input HDF5 file and collects all necessary meta data
int open_hdf5(char *filename, header* h, dedisp_float*& input, int verbose);

// Usage help
void usage(char *prg);
// parsing command line arguments
int parse_cmdline(int argc, char *argv[], cmdline* cmd);
