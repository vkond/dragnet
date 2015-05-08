#include "dragnet_gpu.h"

// Read SIGPROC filterbank header
int read_sigproc_header(FILE *file, header* h) {
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

// Get telescope name by its id
static char *telescope_name(int telescope_id) {

   char *telescope, string[80];
   switch (telescope_id) {
   case 0:
      strcpy(string, "Fake");
      break;
   case 1:
      strcpy(string, "Arecibo");
      break;
   case 2:
      strcpy(string, "Ooty");
      break;
   case 3:
      strcpy(string, "Nancay");
      break;
   case 4:
      strcpy(string, "Parkes");
      break;
   case 5:
      strcpy(string, "Jodrell");
      break;
   case 6:
      strcpy(string, "GBT");
      break;
   case 7:
      strcpy(string, "GMRT");
      break;
   case 8:
      strcpy(string, "Effelsberg");
      break;
   case 9:
      strcpy(string, "ATA");
      break;
   case 10:
      strcpy(string, "UTR-2");
      break;
   case 11:
      strcpy(string, "LOFAR");
      break;
   case 12:
      strcpy(string, "FR606");
      break;
   case 13:
      strcpy(string, "DE601");
      break;
   case 14:
      strcpy(string, "UK608");
      break;
   default:
      strcpy(string, "???????");
      break;
   }
   telescope = (char *) calloc(strlen(string) + 1, 1);
   strcpy(telescope, string);
   return telescope;
}

// Get backend name by its id
static char *backend_name(int machine_id) {

   char *backend, string[80];
   switch (machine_id) {
   case 0:
      strcpy(string, "FAKE");
      break;
   case 1:
      strcpy(string, "PSPM");
      break;
   case 2:
      strcpy(string, "WAPP");
      break;
   case 3:
      strcpy(string, "AOFTM");
      break;
   case 4:
      strcpy(string, "BPP");
      break;
   case 5:
      strcpy(string, "OOTY");
      break;
   case 6:
      strcpy(string, "SCAMP");
      break;
   case 7:
      strcpy(string, "SPIGOT");
      break;
   case 10:
      strcpy(string, "ARTEMIS");
      break;
   case 11:
      strcpy(string, "Cobalt");
      break;
   default:
      strcpy(string, "????");
      break;
   }
   backend = (char *) calloc(strlen(string) + 1, 1);
   strcpy(backend, string);
   return backend;
}
 
// Writes out .inf file
void writeinf(header *h, char *outstem, float dm) {

        char outname[1024];
        char tmp1[100], tmp2[100];
        int itmp;
        int ra_h, ra_m, dec_d, dec_m;
        double ra_s, dec_s;
        FILE *infofile;

        sprintf(outname, "%s_DM%.3f.inf", outstem, dm);

        // first check if file already exists. If it does, then return
        // struct stat info;
        // if (stat(outname, &info) == 0) return;

        if ((infofile=fopen(outname, "w")) == NULL) {
                fprintf(stderr, "Error opening output inf-file!\n");
                exit(1);
        }

        fprintf(infofile, " Data file name without suffix          =  %s_DM%.3f\n", outstem, dm);
        fprintf(infofile, " Telescope used                         =  %s\n", telescope_name(h->telescope_id));
        fprintf(infofile, " Instrument used                        =  %s\n", backend_name(h->machine_id));
        fprintf(infofile, " Object being observed                  =  %s\n", h->source_name);
        ra_h = (int) floor(h->src_raj / 10000.0);
        ra_m = (int) floor((h->src_raj - ra_h * 10000) / 100.0);
        ra_s = h->src_raj - ra_h * 10000 - ra_m * 100;
        dec_d = (int) floor(fabs(h->src_dej) / 10000.0);
        dec_m = (int) floor((fabs(h->src_dej) - dec_d * 10000) / 100.0);
        dec_s = fabs(h->src_dej) - dec_d * 10000 - dec_m * 100;
        if (h->src_dej < 0.0) dec_d = -dec_d;
        fprintf(infofile, " J2000 Right Ascension (hh:mm:ss.ssss)  =  %02d:%02d:%02f\n", ra_h, ra_m, ra_s);
        fprintf(infofile, " J2000 Declination     (dd:mm:ss.ssss)  =  %02d:%02d:%s%f\n", dec_d, dec_m, dec_s < 10 ? "0" : "", dec_s);
        fprintf(infofile, " Data observed by                       =  Unknown\n");
        sprintf(tmp1, "%.15f", h->tstart - (int) floor(h->tstart));
        sscanf(tmp1, "%d.%s", &itmp, tmp2);
	fprintf(infofile, " Epoch of observation (MJD)             =  %d.%s\n", (int) floor(h->tstart), tmp2);
	fprintf(infofile, " Barycentered?           (1=yes, 0=no)  =  0\n");
	fprintf(infofile, " Number of bins in the time series      =  %d\n", h->nsamp);
	fprintf(infofile, " Width of each time series bin (sec)    =  %.15g\n", h->tsamp);
	fprintf(infofile, " Any breaks in the data? (1=yes, 0=no)  =  0\n");
	fprintf(infofile, " Type of observation (EM band)          =  Radio\n");
	fprintf(infofile, " Beam diameter (arcsec)                 =  3600\n");
	fprintf(infofile, " Dispersion measure (cm-3 pc)           =  %.12g\n", dm);
	fprintf(infofile, " Central freq of low channel (Mhz)      =  %.12g\n", h->fch1 - (h->nchan - 1) * fabs(h->foff));
	fprintf(infofile, " Total bandwidth (Mhz)                  =  %.12g\n", fabs(h->foff) * h->nchan);
	fprintf(infofile, " Number of channels                     =  %d\n", h->nchan);
	fprintf(infofile, " Channel bandwidth (Mhz)                =  %.12g\n", fabs(h->foff));
	fprintf(infofile, " Data analyzed by                       =  Unknown\n");
	fprintf(infofile, " Any additional notes:\n    Input filterbank samples have %d bits.\n", h->nbit);

        fclose(infofile);
}
