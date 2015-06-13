#include "dragnet.h"
#include "lofarhdf5.h"
#include "sigproc.h"
#include "mask/mask.h"
#include "mask/vectors.h"

/* ---------- m a i n -----------*/ 

int main(int argc,char *argv[])
{
  int out, next;
  unsigned int i;
  header h;
  cmdline opts;
  void *input = 0;
  void *raw;  // file descriptor of input raw data
  dedisp_bool *killmask = NULL;
  dedisp_float *output=0, *finput = 0;
  dedisp_float *dmlist;
  dedisp_plan plan;
  dedisp_size dm_count, max_delay;
  int64_t nsamp_computed;
  dedisp_error error;
  dedisp_size nbits = 32;
  char *filename, outname[1024];
  // rfi masking
  mask obsmask;  // PRESTO type mask struct
  unsigned int is_mask_apply = 0;  // if we want to apply the mask
  float *padvals = NULL;  // padding values for where RFI is (80% of average)
  int *maskchans = NULL;  // masked channels
  int numzapchan, *zapchan = NULL; // user-defined list of channels to zap

  // initializing cmdline structure
  opts.device_id = 0;
  opts.verbose = 1;
  opts.blocksize = 0;
  strcpy(opts.prefix, "test");
  strcpy(opts.format, "sigproc");
  strcpy(opts.maskfile, "\0");
  strcpy(opts.zapchan, "\0");
  opts.dm_start = 0.0;
  opts.dm_end = 50.0;
  opts.pulse_width = 4.0;
  opts.dm_tol = 1.25;
  opts.clip_sigma = 0.0;

  // parsing cmdline
  if (argc == 1) usage(argv[0]);
  next = parse_cmdline(argc, argv, &opts);
  if (argc < next + 1) { 
   fprintf(stderr, "ERROR: No input file!\n");
   return -1;
  }
  filename = argv[next];

  // open the file
  if (raw_open(filename, opts.format, &h, opts.verbose, raw) != 0) 
   exit(-1);

  /* Get list of user-zapped channels */
  if (strcmp(opts.zapchan, "\0") != 0) {
     zapchan = ranges_to_ivect(opts.zapchan, 0, h.nchan - 1, &numzapchan);
     printf ("Number of user-defined channels to zap: %d (%s)\n", numzapchan, opts.zapchan);
  }

  /* Read an input mask if wanted */
  if (strcmp(opts.maskfile, "\0") != 0) {
    is_mask_apply = 1;
    maskchans = gen_ivect(h.nchan);
    read_mask(opts.maskfile, &obsmask);
    printf("Read mask information from '%s'\n\n", opts.maskfile);
    padvals = (float *)malloc(h.nchan * sizeof(float));
    memset(padvals, 0, h.nchan * sizeof(float));
    determine_padvals(opts.maskfile, &obsmask, padvals);
    // inverse channels in the mask
    inverse_mask(&obsmask, padvals);

    if (h.nbit / 8 == 1) {
      finput=(dedisp_float *)malloc((opts.blocksize > h.nsamp ? h.nsamp : opts.blocksize) * h.nchan * sizeof(dedisp_float));
      if (finput==NULL) {
         fprintf(stderr, "ERROR: Failed to allocate finput array\n");
         return -1;
      }
    }
   } else { obsmask.numchan = obsmask.numint = 0; }

  // checking the block size
  if (opts.blocksize <= 0 || opts.blocksize > h.nsamp) opts.blocksize = h.nsamp;

  // Intialize GPU
  if (opts.verbose) printf("\nIntializing GPU (device %d)\n", opts.device_id);
  error = dedisp_set_device(opts.device_id);
  if (error != DEDISP_NO_ERROR) {
    printf("ERROR: Could not set GPU device: %s\n", dedisp_get_error_string(error));
    return -1;
  }

  // Create a dedispersion plan
  if (opts.verbose) printf("Creating dedispersion plan...\n");
  error = dedisp_create_plan(&plan, h.nchan, h.tsamp, h.fch1, h.foff);
  if (error != DEDISP_NO_ERROR) {
    printf("ERROR: Could not create dedispersion plan: %s\n", dedisp_get_error_string(error));
    return -1;
  }

  // Generate a list of dispersion measures for the plan
  if (opts.dm_step == 0) {
    if (opts.verbose) printf("Generating optimal DM trials...\n");
    error = dedisp_generate_dm_list(plan, opts.dm_start, opts.dm_end, opts.pulse_width, opts.dm_tol);
    if (error != DEDISP_NO_ERROR) {
      printf("ERROR: Failed to generate DM list: %s\n", dedisp_get_error_string(error));
      return -1;
    }
  } else {
    if (opts.verbose) printf("Generating linear DM trials...\n");
    dm_count=(int) ceil((opts.dm_end - opts.dm_start)/opts.dm_step) + 1;
    dmlist=(dedisp_float *) calloc(sizeof(dedisp_float), dm_count);
    for (i=0; i<dm_count; i++) {
      dmlist[i]=(dedisp_float) opts.dm_start + opts.dm_step*i;
//      printf("dm[%d] = %f\n", i, dmlist[i]);
    }
    error=dedisp_set_dm_list(plan, dmlist, dm_count);
    if (error != DEDISP_NO_ERROR) {
      printf("ERROR: Failed to generate DM list: %s\n", dedisp_get_error_string(error));
      return -1;
    }
  }

  // zapping given channels using Dedisp's killmask
  // freq order is: first chan is high freq
  // VLAD: for whatever reason I couldn't make it work
  // I tried different freqs order and inversed boolean values (0|1)
  /*
  if (strcmp(opts.zapchan, "\0") != 0) {
    killmask = (dedisp_bool *)malloc(sizeof(dedisp_bool) * h.nchan);
    memset(killmask, 1, sizeof(dedisp_bool) * h.nchan);
    for (int jj = 0; jj < numzapchan; jj++) killmask[h.nchan - 1 - zapchan[jj]] = 0;
    error = dedisp_set_killmask(plan, killmask);
    if (error != DEDISP_NO_ERROR) {
      printf("ERROR: Failed to set killmask: %s\n", dedisp_get_error_string(error));
      return -1;
    }
  }
  */

  // Get specifics of the computed dedispersion plan
  dm_count = dedisp_get_dm_count(plan);
  max_delay = dedisp_get_max_delay(plan);
  nsamp_computed = h.nsamp-max_delay;
  dmlist=(dedisp_float *)dedisp_get_dm_list(plan);

  // checking if our blocksize is smaller than max_delay
  if (opts.blocksize <= max_delay) {
      printf("ERROR: input data blocksize (%lld) is smaller than Max DM delay (%ld)!\n", opts.blocksize, max_delay);
      return -1;
  }

  // Print information
  if (opts.verbose) {
    printf("----------------------------- DM COMPUTATIONS  ----------------------------\n");
    printf("Computing %ld DMs from %f to %f pc/cc\n", dm_count, dmlist[0], dmlist[dm_count-1]);
    printf("Max DM delay is %ld samples (%f seconds)\n", max_delay, max_delay*h.tsamp);
    printf("Computing %ld out of %d total samples (%.2f%% efficiency)\n", nsamp_computed, h.nsamp, 100.0*(dedisp_float)nsamp_computed/h.nsamp);
    if (opts.dm_step==0.0) printf("Pulse width: %f, DM tolerance: %f\n", opts.pulse_width, opts.dm_tol);
    printf("Output data array size : %ld MB\n", (dm_count * nsamp_computed * (nbits/8))/(1<<20));
    printf("\n");
  }

  // Setting maximum gulp_size
  dedisp_size gulp_max, gulp_size = dedisp_get_gulp_size(plan);
  if (gulp_size < opts.blocksize - max_delay) gulp_max = gulp_size;
  else gulp_max = opts.blocksize - max_delay;
  printf("Setting gulp_size from %d to %d\n", (int)gulp_size, (int)gulp_max);
  error = dedisp_set_gulp_size(plan, gulp_max);
  //printf("Current gulp_size = %d\n", (int)dedisp_get_gulp_size(plan));
  if (error != DEDISP_NO_ERROR) {
    printf("ERROR: Failed to set gulp_size: %s\n", dedisp_get_error_string(error));
    return -1;
  }

  // Loop over data blocks

  int idata = 0; // loop counter
  int64_t read_samples, to_read, isamp_computed = 0;
  clock_t startclock;

  do {
      if (opts.verbose) printf("Data block: %d\n", idata);

      to_read = (isamp_computed + opts.blocksize > h.nsamp ? h.nsamp - isamp_computed : opts.blocksize);
      // Reading the input data
      read_samples = raw_read(to_read, max_delay, &h, input, raw);
   
      // to zap given channels
      // VLAD: more effectively it would be to use Dedisp's killmask, 
      // but I couldn't make it work (see above)
      if (strcmp(opts.zapchan, "\0") != 0) {
        dedisp_byte *ptr = (dedisp_byte *)input;
        for (int64_t ii=0; ii<to_read; ii++)
          for (int jj=0; jj<numzapchan; jj++) ptr[ii * h.nchan + h.nchan - 1 - zapchan[jj]] = 0;
      }

      // applying rfi mask
      if (is_mask_apply) {
          // filling in finput array
          if (h.nbit / 8 == 1) {
             dedisp_byte *ptr = (dedisp_byte *)input;
             // I need to re-arrange the order of channels
             // as Presto assumed lowest channel to be the first
             for (int64_t ii = 0; ii<to_read * h.nchan; ii++) finput[ii] = ptr[ii];
          } else finput = (dedisp_float *)input;
          apply_mask(finput, &h, to_read, isamp_computed, opts.clip_sigma, padvals, maskchans, &obsmask);
      }

      // Allocate space for the output data
      output=(dedisp_float *)realloc(output, (to_read - max_delay) * dm_count * nbits/8);
      if (output==NULL) {
        fprintf(stderr, "ERROR: Failed to allocate output array\n");
        return -1;
      }

      // Perform computation
      if (opts.verbose) printf("Dedispersing on the GPU\n");
      startclock=clock();
      error = dedisp_execute(plan, to_read, is_mask_apply == 1 ? (dedisp_byte *)finput : (dedisp_byte *)input, \
              is_mask_apply == 1 ? 32 : h.nbit, (dedisp_byte *)output, nbits, DEDISP_USE_DEFAULT);
      if (error != DEDISP_NO_ERROR) {
        fprintf(stderr, "ERROR: Failed to execute dedispersion plan: %s\n", dedisp_get_error_string(error));
        return -1;
      }
      if (opts.verbose) printf("Dedispersion took %.2f seconds\n",(double)(clock()-startclock)/CLOCKS_PER_SEC);

      // Write output DM trials
      for (i=0; i<dm_count; i++) {
        // Generate output file name
        sprintf(outname, "%s_DM%.3f.dat", opts.prefix, dmlist[i]);
        if ( (out = open(outname, (idata == 0 ? O_TRUNC : O_APPEND) | O_CREAT | O_WRONLY | O_LARGEFILE, 0644)) == -1 ) {
          fprintf(stderr, "Error opening %s\n", outname);
          return -1;
        }

        // Write buffer
        write(out, output + i*(to_read - max_delay), sizeof(dedisp_float)*(to_read - max_delay));
    
        // Close file
        close(out);

        // Write inf file
        if (idata == 0) writeinf(&h, opts.prefix, dmlist[i], max_delay);
      }

      isamp_computed += read_samples;
      idata++;
  } while (isamp_computed < nsamp_computed);

  if (opts.verbose) printf("\nFinish.\n");

  // Clean up
  if (is_mask_apply) {
    free(padvals);
    free_mask(obsmask);
    vect_free(maskchans);
    if (h.nbit / 8 == 1) {
     if (finput != NULL) free(finput);
    }
  }
  if (killmask != NULL) free(killmask);
  if (zapchan != NULL) free(zapchan);
  if (input != NULL) free(input);
  if (output != NULL) free(output);
  dedisp_destroy_plan(plan);
  raw_close(raw);

  return 0;
}
