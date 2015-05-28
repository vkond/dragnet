#include "sigproc.h"

// constructor
Sigproc::Sigproc(char *filename, header* h, int verbose) : Input() {
  // reading Sigproc input file
  if ((open(filename, h, verbose)) != 0) {
     fprintf(stderr, "Error reading file %s\n", filename);
     exit(-1);
  }
}

// closing the file
void Sigproc::close() {
    fclose(fd);
    if (input != NULL) free(input);
}

// destructor
Sigproc::~Sigproc() {
    close();
}

// opens input Sigproc file and read the header
int Sigproc::open(char *filename, header* h, int verbose) {

  current_sample = 0;
  
  // reading header and data from the SIGPROC file
  if ((fd = fopen(filename, "rb")) == NULL) return -1;
  // reading Sigproc header
  if ((read_sigproc_header(fd, h)) != 0) {
    fprintf(stderr, "Error reading header from file %s\n", filename);
    return -1;
  }
  // Print information
  if (verbose) {
    printf("----------------------------- INPUT DATA ---------------------------------\n");
    printf("Frequency of highest channel              : %f MHz\n", h->fch1);
    printf("Bandwidth                                 : %f MHz\n", fabs(h->foff)*h->nchan);
    printf("Number of channels (channel width)        : %d (%f MHz)\n", h->nchan, fabs(h->foff));
    printf("Sample time                               : %f us\n", h->tsamp*1e6);
    printf("Observation duration                      : %f s (%d samples)\n", h->tsamp*h->nsamp, h->nsamp);
    printf("Number of polarizations/bit depth         : %d/%d\n", h->nif, h->nbit);
    printf("Input data array size                     : %lu MB\n", h->buffersize/(1<<20));
    printf("Header size                               : %lu bytes\n", h->headersize);
    printf("\n");
  }

  // Exit on wrong type of input data
  if (h->nif != 1) {
    fprintf(stderr, "Wrong number of polarizations (not 1). Exiting.\n");
    return -1;
  }
  if (h->nbit != 8) {
    fprintf(stderr, "Wrong bit depth (not 8). Exiting.\n");
    return -1;
  }
  return 0;
}

// Read SIGPROC filterbank header
int Sigproc::read_sigproc_header(FILE *file, header* h) {
  int nchar, nbytes=0, expecting_fchannel = 0;
  char string[80];
  float temp;

  // Read header parameters
  for (;;) {
    // Read string size
    strcpy(string,"ERROR");
    fread(&nchar, sizeof(int), 1, file);

    // Skip wrong strings
    if (!(nchar > 1 && nchar < 80)) continue;

    // Increate byte counter
    nbytes += nchar;

    // Read string
    fread(string, nchar, 1, file);
    string[nchar]='\0';

    // Exit at end of header
    if (strcmp(string, "HEADER_END") == 0) break;

    // Read parameters
    if (strcmp(string, "tsamp") == 0) fread(&(h->tsamp), sizeof(double), 1, file);
    else if (strcmp(string,"tstart") == 0) fread(&(h->tstart), sizeof(double), 1, file);
    else if (strcmp(string,"fch1") == 0) fread(&(h->fch1), sizeof(double), 1, file);
    else if (strcmp(string,"foff") == 0) fread(&(h->foff), sizeof(double), 1, file);
    else if (strcmp(string,"nchans") == 0) fread(&(h->nchan), sizeof(int), 1, file);
    else if (strcmp(string,"nifs") == 0) fread(&(h->nif), sizeof(int), 1, file);
    else if (strcmp(string,"nbits") == 0) fread(&(h->nbit), sizeof(int), 1, file);
    else if (strcmp(string,"nsamples") == 0) fread(&(h->nsamp), sizeof(int), 1, file);
    else if (strcmp(string,"az_start") == 0) fread(&(h->az_start), sizeof(double), 1, file);
    else if (strcmp(string,"za_start") == 0) fread(&(h->za_start), sizeof(double), 1, file);
    else if (strcmp(string,"src_raj") == 0) fread(&(h->src_raj), sizeof(double), 1, file);
    else if (strcmp(string,"src_dej") == 0) fread(&(h->src_dej), sizeof(double), 1, file);
    else if (strcmp(string,"telescope_id") == 0) fread(&(h->telescope_id), sizeof(int), 1, file);
    else if (strcmp(string,"machine_id") == 0) fread(&(h->machine_id), sizeof(int), 1, file);
    else if (strcmp(string,"nbeams") == 0) fread(&(h->nbeam), sizeof(int), 1, file);
    else if (strcmp(string,"ibeam") == 0) fread(&(h->ibeam), sizeof(int), 1, file);
    else if (strcmp(string,"source_name") == 0) strcpy(h->source_name, string);
    else if (strcmp(string, "rawdatafile") == 0) strcpy(h->inpfile, string);
    else if (strcmp(string, "FREQUENCY_START") == 0) expecting_fchannel = 1;
    else if (strcmp(string, "FREQUENCY_END") == 0) expecting_fchannel = 0;
    else if ((strcmp(string, "fchannel") == 0) && expecting_fchannel) {
            if (expecting_fchannel == 3) fread(&temp, sizeof(double), 1, file);
            if (expecting_fchannel == 2) {
              fread(&(h->foff), sizeof(double), 1, file);
              h->foff = h->foff - h->fch1;
              expecting_fchannel = 3;
            }
            if (expecting_fchannel == 1) {
              fread(&(h->fch1), sizeof(double), 1, file);
              expecting_fchannel = 2;
            }
    }
  }

  // Get header and buffer sizes
  h->headersize = (int64_t)ftell(file);
  fseek(file, 0, SEEK_END);
  h->buffersize = ftell(file) - h->headersize;
  h->nsamp = h->buffersize/(h->nchan * h->nif * h->nbit/8);

  // Reset file pointer to start of buffer
  rewind(file);
  fseek(file, h->headersize, SEEK_SET);

  return 0;
}

// reading the nsamples from the file
int64_t Sigproc::read(int64_t nsamples, int64_t shift_back, header* h, void*& out) {

  // how many samples we actually read
  int64_t read_samples = current_sample + nsamples > h->nsamp ? h->nsamp - current_sample : nsamples;
  // allocating memory for input buffer
  input=(dedisp_byte *) realloc(input, sizeof(dedisp_byte) * read_samples * h->nchan);

  fread(input, sizeof(dedisp_byte), read_samples * h->nchan, fd);

  out = input; // set the (void*) pointer to where the data are
  current_sample += (read_samples - shift_back);
  // rewinding back by shift_back samples
  fseek (fd, -shift_back * (h->nbit/8) * h->nif * h->nchan, SEEK_CUR);
  return read_samples - shift_back;
}
