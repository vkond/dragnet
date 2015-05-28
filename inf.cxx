#include "dragnet.h"

// Get telescope name by its id
char *telescope_name(int telescope_id) {

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
char *backend_name(int machine_id) {

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
void writeinf(header *h, char *outstem, float dm, int64_t shift_back) {

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
	fprintf(infofile, " Number of bins in the time series      =  %d\n", h->nsamp - (int)shift_back);
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
