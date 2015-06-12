#include "dragnet.h"
#include "input.h"
#include "sigproc.h"
#include "lofarhdf5.h"

// open input file
int raw_open(char *filename, char *format, header *h, int verbose, void*& fd) {

  // checking the data format
  if ((strcmp(format, "sigproc")) == 0) {
   // opening the input Sigproc file
   fd = new Sigproc(filename, h, verbose);
  } else if ((strcmp(format, "hdf5")) == 0) {
   // opening the input HDF5 file
   fd = new HDF5(filename, h, verbose);
  } else {
   fprintf(stderr, "ERROR: Can't recognise the format of input data!\n");
   return -1;
  }
 return 0;
}

// read input file
int64_t raw_read(int64_t nsamples, int64_t shift_back, header* h, void*& out, void*& fd) {
 int64_t read_samples = ((Input *)fd)->read(nsamples, shift_back, h, out);
 return read_samples;
}

// close input file
void raw_close(void*& fd) {
 ((Input *)fd)->close();
}
