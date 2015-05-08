#include "dragnet_gpu.h"

/* ---------- m a i n -----------*/ 

int main(int argc,char *argv[])
{
  int i, out, next;
  header h;
  cmdline opts;
  void *input = 0;
  dedisp_float *output=0;
  dedisp_float *dmlist;
  dedisp_plan plan;
  dedisp_size dm_count, max_delay, nsamp_computed;
  dedisp_error error;
  dedisp_size nbits = 32;
  char *filename, outname[1024];

  // initializing cmdline structure
  opts.device_id = 0;
  opts.verbose = 1;
  strcpy(opts.prefix, "test");
  strcpy(opts.format, "sigproc");
  opts.dm_start = 0.0;
  opts.dm_end = 50.0;
  opts.dm_start = 0.0;
  opts.pulse_width = 4.0;
  opts.dm_tol = 1.25;

  // parsing cmdline
  if (argc == 1) usage(argv[0]);
  next = parse_cmdline(argc, argv, &opts);
  if (argc < next + 1) { 
   fprintf(stderr, "ERROR: No input file!\n");
   return -1;
  }
  filename = argv[next];


  // checking the data format
  if ((strcmp(opts.format, "sigproc")) == 0) {

    dedisp_byte *tmpinput=0;
    FILE *file;

    // reading header and data from the SIGPROC file
    if ((file = fopen(filename, "rb")) == NULL) {
     fprintf(stderr, "Error reading file %s\n", filename);
     return -1;
    }
    // reading Sigproc header
    if ((read_sigproc_header(file, &h)) != 0) {
     fprintf(stderr, "Error reading header from file %s\n", filename);
     return -1;
    }

    // Print information
    if (opts.verbose) {
      printf("----------------------------- INPUT DATA ---------------------------------\n");
      printf("Frequency of highest channel              : %f MHz\n", h.fch1);
      printf("Bandwidth                                 : %f MHz\n", fabs(h.foff)*h.nchan);
      printf("Number of channels (channel width)        : %d (%f MHz)\n", h.nchan,fabs(h.foff));
      printf("Sample time                               : %f us\n", h.tsamp*1e6);
      printf("Observation duration                      : %f s (%d samples)\n", h.tsamp*h.nsamp,h.nsamp);
      printf("Number of polarizations/bit depth         : %d/%d\n", h.nif, h.nbit);
      printf("Input data array size                     : %lu MB\n", h.buffersize/(1<<20));
      printf("Header size                               : %lu bytes\n", h.headersize);
      printf("\n");
    }

    // Exit on wrong type of input data
    if (h.nif != 1) {
      fprintf(stderr, "Wrong number of polarizations (not 1). Exiting.\n");
      return -1;
    }
    if (h.nbit != 8) {
      fprintf(stderr, "Wrong bit depth (not 8). Exiting.\n");
      return -1;
    }

    // Read buffer
    if (opts.verbose) printf("Reading file\n");
    tmpinput=(dedisp_byte *) malloc(sizeof(dedisp_byte)*h.buffersize);
    fread(tmpinput, sizeof(dedisp_byte), h.buffersize, file);
    input = tmpinput;

    // Close file;
    fclose(file);

  } 
  else if ((strcmp(opts.format, "hdf5")) == 0) {

    dedisp_float *tmpinput=0;

    // reading HDF5 input file
    if ((open_hdf5(filename, &h, tmpinput, opts.verbose)) != 0) {
     fprintf(stderr, "Error reading file %s\n", filename);
     return -1;
    }
    input = tmpinput;
  }
  else {
   fprintf(stderr, "ERROR: Can't recognise the format of input data!\n");
   return -1;
  }




  // Intialize GPU
  if (opts.verbose) printf("\nIntializing GPU (device %d)\n", opts.device_id);
  error = dedisp_set_device(opts.device_id);
  if (error != DEDISP_NO_ERROR) {
    printf("ERROR: Could not set GPU device: %s\n", dedisp_get_error_string(error));
    return -1;
  }

  // Create a dedispersion plan
  if (opts.verbose) printf("Creating dedispersion plan\n");
  error = dedisp_create_plan(&plan, h.nchan, h.tsamp, h.fch1, h.foff);
  if (error != DEDISP_NO_ERROR) {
    printf("ERROR: Could not create dedispersion plan: %s\n", dedisp_get_error_string(error));
    return -1;
  }

  // Generate a list of dispersion measures for the plan
  if (opts.dm_step == 0) {
    if (opts.verbose) printf("Generating optimal DM trials\n");
    error = dedisp_generate_dm_list(plan, opts.dm_start, opts.dm_end, opts.pulse_width, opts.dm_tol);
    if (error != DEDISP_NO_ERROR) {
      printf("ERROR: Failed to generate DM list: %s\n", dedisp_get_error_string(error));
      return -1;
    }
  } else {
    // Generate a list of dispersion measures for the plan
    if (opts.verbose) printf("Generating linear DM trials\n");
    dm_count=(int) ceil((opts.dm_end - opts.dm_start)/opts.dm_step) + 1;
    dmlist=(dedisp_float *) calloc(sizeof(dedisp_float), dm_count);
    for (i=0; i<dm_count; i++) {
      dmlist[i]=(dedisp_float) opts.dm_start + opts.dm_step*i;
      printf("dm[%d] = %f\n", i, dmlist[i]);
    }
    error=dedisp_set_dm_list(plan, dmlist, dm_count);
    if (error != DEDISP_NO_ERROR) {
      printf("ERROR: Failed to generate DM list: %s\n", dedisp_get_error_string(error));
      return -1;
    }
  }

  // Get specifics of the computed dedispersion plan
  dm_count = dedisp_get_dm_count(plan);
  max_delay = dedisp_get_max_delay(plan);
  nsamp_computed = h.nsamp-max_delay;
  dmlist=(dedisp_float *)dedisp_get_dm_list(plan);

  // Print information
  if (opts.verbose) {
    printf("----------------------------- DM COMPUTATIONS  ----------------------------\n");
    printf("Computing %ld DMs from %f to %f pc/cm^3\n", dm_count, dmlist[0], dmlist[dm_count-1]);
    printf("Max DM delay is %ld samples (%f seconds)\n", max_delay, max_delay*h.tsamp);
    printf("Computing %ld out of %d total samples (%.2f%% efficiency)\n", nsamp_computed, h.nsamp, 100.0*(dedisp_float)nsamp_computed/h.nsamp);
    if (opts.dm_step==0.0) printf("Pulse width: %f, DM tolerance: %f\n", opts.pulse_width, opts.dm_tol);
    printf("Output data array size : %ld MB\n", (dm_count * nsamp_computed * (nbits/8))/(1<<20));
    printf("\n");
  }

  // Allocate space for the output data
  output=(dedisp_float *)malloc(nsamp_computed * dm_count * nbits/8);
  if (output==NULL) {
    fprintf(stderr, "ERROR: Failed to allocate output array\n");
    return -1;
  }

  // Perform computation
  if (opts.verbose) printf("Dedispersing on the GPU\n");
  clock_t startclock=clock();
  error = dedisp_execute(plan, h.nsamp, (dedisp_byte *)input, h.nbit, (dedisp_byte *)output, nbits, DEDISP_USE_DEFAULT);
  if (error != DEDISP_NO_ERROR) {
    fprintf(stderr, "ERROR: Failed to execute dedispersion plan: %s\n", dedisp_get_error_string(error));
    return -1;
  }
  if (opts.verbose) printf("Dedispersion took %.2f seconds\n",(double)(clock()-startclock)/CLOCKS_PER_SEC);

  // Write output DM trials
  for (i=0; i<dm_count; i++) {
    // Generate output file name
    sprintf(outname, "%s_DM%.3f.dat", opts.prefix, dmlist[i]);
    if ( (out = open(outname, O_TRUNC | O_CREAT | O_WRONLY | O_LARGEFILE, 0644)) == -1 ) {
      fprintf(stderr, "Error opening %s\n", outname);
      return -1;
    }

    // Write buffer
    write(out, output + i*nsamp_computed, sizeof(dedisp_float)*nsamp_computed);
    
    // Close file
    close(out);

    // Write inf file
    writeinf(&h, opts.prefix, dmlist[i]);
  }

  // Clean up
  if (input != NULL) free(input);
  free(output);
  dedisp_destroy_plan(plan);

  return 0;
}
