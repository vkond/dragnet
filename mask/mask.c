#include <string.h>
#include <math.h>
#include "chkio.h"
#include "vectors.h"
#include "mask.h"
#include "../dragnet.h"

#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

static int find_num(int num, int *arr, int arrlen);
static int merge_no_dupes(int *arr1, int len1, int *arr2, int len2, int *merged);

/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipies in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 */

/* Fast computation of the median of an array. */
/* Note:  It messes up the order!              */

float median(float arr[], int n)
{
   int low, high;
   int median;
   int middle, ll, hh;

   low = 0;
   high = n - 1;
   median = (low + high) / 2;
   for (;;) {
      if (high <= low)          /* One element only */
         return arr[median];

      if (high == low + 1) {    /* Two elements only */
         if (arr[low] > arr[high])
            ELEM_SWAP(arr[low], arr[high]);
         return arr[median];
      }

      /* Find median of low, middle and high items; swap into position low */
      middle = (low + high) / 2;
      if (arr[middle] > arr[high])
         ELEM_SWAP(arr[middle], arr[high]);
      if (arr[low] > arr[high])
         ELEM_SWAP(arr[low], arr[high]);
      if (arr[middle] > arr[low])
         ELEM_SWAP(arr[middle], arr[low]);

      /* Swap low item (now in position middle) into position (low+1) */
      ELEM_SWAP(arr[middle], arr[low + 1]);

      /* Nibble from each end towards middle, swapping items when stuck */
      ll = low + 1;
      hh = high;
      for (;;) {
         do
            ll++;
         while (arr[low] > arr[ll]);
         do
            hh--;
         while (arr[hh] > arr[low]);

         if (hh < ll)
            break;

         ELEM_SWAP(arr[ll], arr[hh]);
      }

      /* Swap middle item (in position low) back into correct position */
      ELEM_SWAP(arr[low], arr[hh]);

      /* Re-set active partition */
      if (hh <= median)
         low = ll;
      if (hh >= median)
         high = hh - 1;
   }
}

void avg_var(float *x, int n, double *mean, double *var)
/* For a float vector, *x, of length n, this routine  */
/* returns the mean and variance of *x.               */
{
   long i;
   double an = 0.0, an1 = 0.0, dx;

   /*  Modified (29 June 98) C version of the following:        */
   /*  ALGORITHM AS 52  APPL. STATIST. (1972) VOL.21, P.226     */
   /*  Returned values were checked with Mathematica 3.01       */

   if (n < 1) {
      printf("\vVector length must be > 0 in avg_var().  Exiting\n");
      exit(1);
   } else {
      *mean = (double) x[0];
      *var = 0.0;
   }

   for (i = 1; i < n; i++) {
      an = (double) (i + 1);
      an1 = (double) (i);
      dx = (x[i] - *mean) / an;
      *var += an * an1 * dx * dx;
      *mean += dx;
   }

   if (n > 1)
      *var /= an1;

   return;
}

int split_root_suffix(char *input, char **root, char **suffix)
/* This routine splits an input string into a root name */
/* + suffix.  Since it allocates the memory for the     */
/* root and suffix dynamically, the calling program     */
/* must free both "root" and "suffix".                  */
/* If the routine finds a suffix, it returns 1, else 0. */
{
   char *sptr = NULL;
   unsigned int len, rootlen = 0, suffixlen = 0;

   len = strlen(input);
   sptr = strrchr(input, '.');
   if (sptr == NULL) {
      *root = (char *) calloc(len + 1, sizeof(char));
      strncpy(*root, input, len);
      return 0;
   } else {
      rootlen = sptr - input;
      *root = (char *) calloc(rootlen + 1, sizeof(char));
      strncpy(*root, input, rootlen);
      suffixlen = len - rootlen - 1;
      *suffix = (char *) calloc(suffixlen + 1, sizeof(char));
      strncpy(*suffix, sptr + 1, suffixlen);
      return 1;
   }
}

int compare_floats(const void *a, const void *b)
    /* qsort comparison function for floats */
{
       const float *da = (const float *) a;
          const float *db = (const float *) b;

             return (*da > *db) - (*da < *db);
}

int compare_ints(const void *a, const void *b)
    /* qsort comparison function for ints */
{
       const int *da = (const int *) a;
          const int *db = (const int *) b;

             return (*da > *db) - (*da < *db);
}

void fill_mask(double timesigma, double freqsigma, double mjd,
               double dtint, double lofreq, double dfreq,
               int numchan, int numint, int ptsperint,
               int num_zap_chans, int *zap_chans, int num_zap_ints,
               int *zap_ints, unsigned char **bytemask, mask * obsmask)
/* Fill a mask structure with the appropriate values */
{
   int ii, jj, count;

   obsmask->timesigma = timesigma;
   obsmask->freqsigma = freqsigma;
   obsmask->mjd = mjd;
   obsmask->dtint = dtint;
   obsmask->lofreq = lofreq;
   obsmask->dfreq = dfreq;
   obsmask->numchan = numchan;
   obsmask->numint = numint;
   obsmask->ptsperint = ptsperint;
   obsmask->num_zap_chans = num_zap_chans;
   if (obsmask->num_zap_chans) {
      obsmask->zap_chans = gen_ivect(obsmask->num_zap_chans);
      for (ii = 0; ii < obsmask->num_zap_chans; ii++)
         obsmask->zap_chans[ii] = zap_chans[ii];
   }
   obsmask->num_zap_ints = num_zap_ints;
   if (obsmask->num_zap_ints) {
      obsmask->zap_ints = gen_ivect(obsmask->num_zap_ints);
      for (ii = 0; ii < obsmask->num_zap_ints; ii++)
         obsmask->zap_ints[ii] = zap_ints[ii];
   }
   obsmask->num_chans_per_int = gen_ivect(obsmask->numint);
   obsmask->chans = (int **) malloc(obsmask->numint * sizeof(int *));
   for (ii = 0; ii < obsmask->numint; ii++) {
      count = 0;
      /* Count the bad channels first */
      for (jj = 0; jj < obsmask->numchan; jj++)
         if ((bytemask[ii][jj] & BADDATA) | 
             (bytemask[ii][jj] & USERZAP))
            count++;
      obsmask->num_chans_per_int[ii] = count;
      if (count) {
         /* Now determine which channels */
         count = 0;
         obsmask->chans[ii] = gen_ivect(obsmask->num_chans_per_int[ii]);
         for (jj = 0; jj < obsmask->numchan; jj++) {
            if ((bytemask[ii][jj] & BADDATA) | 
                (bytemask[ii][jj] & USERZAP))
               obsmask->chans[ii][count++] = jj;
         }
      }
   }
}


void set_oldmask_bits(mask * oldmask, unsigned char **bytemask)
/* Sets the oldmask bit in the appropriate bytes in bytemask */
{
   int ii, jj;

   for (ii = 0; ii < oldmask->numint; ii++)
      for (jj = 0; jj < oldmask->num_chans_per_int[ii]; jj++)
         bytemask[ii][oldmask->chans[ii][jj]] |= OLDMASK;
}


void unset_oldmask_bits(mask * oldmask, unsigned char **bytemask)
/* Unsets the oldmask bits in bytemask */
{
   int ii, jj;

   for (ii = 0; ii < oldmask->numint; ii++)
      for (jj = 0; jj < oldmask->numchan; jj++)
         bytemask[ii][jj] &= ~OLDMASK;
}


void free_mask(mask obsmask)
/* Free the contents of an mask structure */
{
    int ii;

    for (ii = 0; ii < obsmask.numint; ii++) {
        if (obsmask.num_chans_per_int[ii] > 0 &&
            obsmask.num_chans_per_int[ii] <= obsmask.numchan)
            vect_free(obsmask.chans[ii]);
    }
    free(obsmask.chans);
    vect_free(obsmask.num_chans_per_int);
    if (obsmask.num_zap_chans)
        vect_free(obsmask.zap_chans);
    if (obsmask.num_zap_ints)
        vect_free(obsmask.zap_ints);
}


void read_mask(char *maskfilenm, mask * obsmask)
/* Read the contents of a mask structure from a file */
{
   FILE *infile;
   int ii;

   infile = chkfopen(maskfilenm, "rb");
   chkfread(&(obsmask->timesigma), sizeof(double), 1, infile);
   chkfread(&(obsmask->freqsigma), sizeof(double), 1, infile);
   chkfread(&(obsmask->mjd), sizeof(double), 1, infile);
   chkfread(&(obsmask->dtint), sizeof(double), 1, infile);
   chkfread(&(obsmask->lofreq), sizeof(double), 1, infile);
   chkfread(&(obsmask->dfreq), sizeof(double), 1, infile);
   chkfread(&(obsmask->numchan), sizeof(int), 1, infile);
   chkfread(&(obsmask->numint), sizeof(int), 1, infile);
   chkfread(&(obsmask->ptsperint), sizeof(int), 1, infile);
   chkfread(&(obsmask->num_zap_chans), sizeof(int), 1, infile);
   if (obsmask->num_zap_chans) {
      obsmask->zap_chans = gen_ivect(obsmask->num_zap_chans);
      chkfread(obsmask->zap_chans, sizeof(int), obsmask->num_zap_chans, infile);
   }
   chkfread(&(obsmask->num_zap_ints), sizeof(int), 1, infile);
   if (obsmask->num_zap_ints) {
      obsmask->zap_ints = gen_ivect(obsmask->num_zap_ints);
      chkfread(obsmask->zap_ints, sizeof(int), obsmask->num_zap_ints, infile);
   }
   obsmask->num_chans_per_int = gen_ivect(obsmask->numint);
   chkfread(obsmask->num_chans_per_int, sizeof(int), obsmask->numint, infile);
   obsmask->chans = (int **) malloc(obsmask->numint * sizeof(int *));
   for (ii = 0; ii < obsmask->numint; ii++) {
      if (obsmask->num_chans_per_int[ii] > 0 &&
          obsmask->num_chans_per_int[ii] < obsmask->numchan) {
         obsmask->chans[ii] = gen_ivect(obsmask->num_chans_per_int[ii]);
         chkfread(obsmask->chans[ii], sizeof(int),
                  obsmask->num_chans_per_int[ii], infile);
      } else if (obsmask->num_chans_per_int[ii] == obsmask->numchan) {
         int jj;
         obsmask->chans[ii] = gen_ivect(obsmask->num_chans_per_int[ii]);
         for (jj = 0; jj < obsmask->numchan; jj++)
            obsmask->chans[ii][jj] = jj;
      }
   }
   fclose(infile);
}


void calc_avgmedstd(float *arr, int numarr, float fraction,
                    int step, float *avg, float *med, float *std)
/* Calculates the median and middle-'fraction' std deviation  */
/* and average of the array 'arr'.  Values are returned in    */
/* 'avg', 'med' and 'std'.  The array is not modified.        */
{
   int ii, jj, len, start;
   float *tmparr;
   double davg, dstd;

   len = (int) (numarr * fraction + 0.5);
   if (len > numarr || len < 0) {
      printf("fraction (%g) out-of-bounds in calc_avgmedstd()\n", fraction);
      exit(1);
   }
   start = (numarr - len) / 2;
   tmparr = gen_fvect(numarr);
   for (ii = 0, jj = 0; ii < numarr; ii++, jj += step)
      tmparr[ii] = arr[jj];
   qsort(tmparr, numarr, sizeof(float), compare_floats);
   avg_var(tmparr + start, len, &davg, &dstd);
   *avg = (float) davg;
   *med = tmparr[numarr / 2];
   *std = sqrt(dstd);
   vect_free(tmparr);
}


int determine_padvals(char *maskfilenm, mask * obsmask, float *padvals)
// Determine reasonable padding values from the rfifind produced
// *.stats file if it is available.  The pre-allocated vector (of
// length numchan) is in padvals.  Return a '1' if the routine used
// the stats file, return 0 if the padding was set to zeros.
{
   FILE *statsfile;
   int ii, numchan, numint, ptsperint, lobin, numbetween;
   float **dataavg, tmp1, tmp2;
   char *statsfilenm, *root, *suffix;

   if (split_root_suffix(maskfilenm, &root, &suffix) == 0) {
      printf("\nThe mask filename (%s) must have a suffix!\n\n", maskfilenm);
      exit(1);
   } else {
      /* Determine the stats file name */
      statsfilenm = (char *)calloc(strlen(maskfilenm) + 2, sizeof(char));
      sprintf(statsfilenm, "%s.stats", root);
      free(root);
      free(suffix);
      /* Check to see if the file exists */
      printf("Attempting to read the data statistics from '%s'...\n", statsfilenm);
      statsfile = chkfopen(statsfilenm, "rb");
      free(statsfilenm);
      if (statsfile) {          /* Read the stats */
         chkfread(&numchan, sizeof(int), 1, statsfile);
         chkfread(&numint, sizeof(int), 1, statsfile);
         chkfread(&ptsperint, sizeof(int), 1, statsfile);
         chkfread(&lobin, sizeof(int), 1, statsfile);
         chkfread(&numbetween, sizeof(int), 1, statsfile);
         dataavg = gen_fmatrix(numint, numchan);
         /* These are the powers */
         chkfread(dataavg[0], sizeof(float), numchan * numint, statsfile);
         /* These are the averages */
         chkfread(dataavg[0], sizeof(float), numchan * numint, statsfile);
         /* Set the padding values equal to the mid-80% channel averages */
         for (ii = 0; ii < numchan; ii++)
            calc_avgmedstd(dataavg[0] + ii, numint, 0.8, numchan,
                           padvals + ii, &tmp1, &tmp2);
         printf
             ("...succeded.  Set the padding values equal to the mid-80%% channel averages.\n");
         vect_free(dataavg[0]);
         vect_free(dataavg);
         fclose(statsfile);
         return 1;
      } else {
          /* This is a temporary solution */
          for (ii = 0; ii < obsmask->numchan; ii++)
              padvals[ii] = 0.0;
          printf("...failed.\n  Set the padding values to 0.\n");
          return 0;
      }
   }
}


void write_mask(char *maskfilenm, mask * obsmask)
/* Write the contents of an mask structure to a file */
{
   FILE *outfile;
   int ii;

   outfile = chkfopen(maskfilenm, "wb");
   chkfwrite(&(obsmask->timesigma), sizeof(double), 1, outfile);
   chkfwrite(&(obsmask->freqsigma), sizeof(double), 1, outfile);
   chkfwrite(&(obsmask->mjd), sizeof(double), 1, outfile);
   chkfwrite(&(obsmask->dtint), sizeof(double), 1, outfile);
   chkfwrite(&(obsmask->lofreq), sizeof(double), 1, outfile);
   chkfwrite(&(obsmask->dfreq), sizeof(double), 1, outfile);
   chkfwrite(&(obsmask->numchan), sizeof(int), 1, outfile);
   chkfwrite(&(obsmask->numint), sizeof(int), 1, outfile);
   chkfwrite(&(obsmask->ptsperint), sizeof(int), 1, outfile);
   chkfwrite(&(obsmask->num_zap_chans), sizeof(int), 1, outfile);
   if (obsmask->num_zap_chans)
      chkfwrite(obsmask->zap_chans, sizeof(int), obsmask->num_zap_chans, outfile);
   chkfwrite(&(obsmask->num_zap_ints), sizeof(int), 1, outfile);
   if (obsmask->num_zap_ints)
      chkfwrite(obsmask->zap_ints, sizeof(int), obsmask->num_zap_ints, outfile);
   chkfwrite(obsmask->num_chans_per_int, sizeof(int), obsmask->numint, outfile);
   for (ii = 0; ii < obsmask->numint; ii++) {
      if (obsmask->num_chans_per_int[ii] > 0 &&
          obsmask->num_chans_per_int[ii] < obsmask->numchan) {
         chkfwrite(obsmask->chans[ii], sizeof(int),
                   obsmask->num_chans_per_int[ii], outfile);
      }
   }
   fclose(outfile);
}


int check_mask(double starttime, double duration, mask * obsmask, int *maskchans)
/* Return value is the number of channels to mask.  The */
/* channel numbers are placed in maskchans (which must  */
/* have a length of numchan).  If -1 is returned, all   */
/* channels should be masked.                           */
{
   int loint, hiint;
   double endtime;
   static int old_loint = -1, old_hiint = -1, old_numchan = 0;

   /*
      static int firsttime = 1;
      if (firsttime){
      int ii;
      printf("\n\n numzapints = %d\n : ", obsmask->num_zap_ints);
      for (ii=0; ii<obsmask->num_zap_ints; ii++)
      printf("%d ", obsmask->zap_ints[ii]);
      printf("\n\n numzapchans = %d\n : ", obsmask->num_zap_chans);
      for (ii=0; ii<obsmask->num_zap_chans; ii++)
      printf("%d ", obsmask->zap_chans[ii]);
      printf("\n\n");
      firsttime = 0;
      }
    */

   endtime = starttime + duration;
   loint = (int) (starttime / obsmask->dtint);
   hiint = (int) (endtime / obsmask->dtint);

   /* Mask the same channels as for the last call */
   if (loint == old_loint && hiint == old_hiint)
       return old_numchan;

   /* Make sure that we aren't past the last interval */
   if (loint >= obsmask->numint)
       loint = obsmask->numint - 1;
   if (hiint >= obsmask->numint)
       hiint = loint;

   if ((loint >= obsmask->numint + 1) ||
       (hiint >= obsmask->numint + 1)) {
       printf("Warning!!  Trying to use a mask interval well after the mask ends!\n");
   }

   /* Determine new channels to mask */
   if (loint == hiint) {
      old_loint = old_hiint = loint;
      /* Check to see if this is an interval where we zap all the channels */
      if (obsmask->num_zap_ints) {
         if (find_num(loint, obsmask->zap_ints, obsmask->num_zap_ints)) {
            old_numchan = -1;
            return old_numchan;
         }
      }
      /* Merge the overall channels to zap with the local channels to zap */
      old_numchan = merge_no_dupes(obsmask->zap_chans,
                                   obsmask->num_zap_chans,
                                   obsmask->chans[loint],
                                   obsmask->num_chans_per_int[loint], maskchans);
   } else {                     /* We are straddling a rfifind interval boundary */
      int *tmpchans;

      old_loint = loint;
      old_hiint = hiint;
      /* Check to see if this is an interval where we zap all the channels */
      if (obsmask->num_zap_ints) {
         if (find_num(loint, obsmask->zap_ints, obsmask->num_zap_ints)) {
            old_numchan = -1;
            return old_numchan;
         }
         if (find_num(hiint, obsmask->zap_ints, obsmask->num_zap_ints)) {
            old_numchan = -1;
            return old_numchan;
         }
      }
      /* Merge the overall channels to zap with the loint channels to zap */
      if (obsmask->num_zap_chans) {
         tmpchans = gen_ivect(obsmask->numchan);
         old_numchan = merge_no_dupes(obsmask->zap_chans,
                                      obsmask->num_zap_chans,
                                      obsmask->chans[loint],
                                      obsmask->num_chans_per_int[loint], tmpchans);
      } else {
         tmpchans = obsmask->zap_chans;
         old_numchan = obsmask->num_zap_chans;
      }
      /* Merge the loint+overall channels to zap with the hiint channels to zap */
      old_numchan = merge_no_dupes(tmpchans,
                                   old_numchan,
                                   obsmask->chans[hiint],
                                   obsmask->num_chans_per_int[hiint], maskchans);
      if (obsmask->num_zap_chans)
         vect_free(tmpchans);
   }
   return old_numchan;
}


static int find_num(int num, int *arr, int arrlen)
{
   int ii;

   /* Note:  I should make sure the array is sorted and do a binary search */
   for (ii = 0; ii < arrlen; ii++)
      if (arr[ii] == num)
         return 1;
   return 0;
}


static int merge_no_dupes(int *arr1, int len1, int *arr2, int len2, int *merged)
{
   int ptr1 = 0, ptr2 = 0, count = 0;

   while (1) {
      if (ptr1 == len1) {
         while (ptr2 < len2)
            merged[count++] = arr2[ptr2++];
         break;
      } else if (ptr2 == len2) {
         while (ptr1 < len1)
            merged[count++] = arr1[ptr1++];
         break;
      }
      if (arr1[ptr1] < arr2[ptr2])
         merged[count++] = arr1[ptr1++];
      else if (arr1[ptr1] > arr2[ptr2])
         merged[count++] = arr2[ptr2++];
      else {
         merged[count++] = arr1[ptr1];
         ptr1++;
         ptr2++;
      }
   }
   return count;
}

/* NEW Clipping Routine (uses channel running averages) */
int clip_times(float *rawdata, int ptsperblk, int numchan,
               float clip_sigma, float *good_chan_levels)
// Perform time-domain clipping of rawdata.  This is a 2D array with
// ptsperblk*numchan points, each of which is a float.  The clipping
// is done at clip_sigma sigma above/below the running mean.  The
// up-to-date running averages of the channels are returned in
// good_chan_levels (which must be pre-allocated).
{
   static float *chan_running_avg;
   static float running_avg = 0.0, running_std = 0.0;
   static int blocksread = 0, firsttime = 1;
   float *zero_dm_block, *ftmp, *powptr;
   double *chan_avg_temp;
   float current_med, trigger;
   double current_avg = 0.0, current_std = 0.0;
   int ii, jj, clipit = 0, clipped = 0;

   if (firsttime) {
       chan_running_avg = gen_fvect(numchan);
       firsttime = 0;
   }
   chan_avg_temp = gen_dvect(numchan);
   zero_dm_block = gen_fvect(ptsperblk);
   ftmp = gen_fvect(ptsperblk);

   /* Calculate the zero DM time series */
   for (ii = 0; ii < ptsperblk; ii++) {
      zero_dm_block[ii] = 0.0;
      powptr = rawdata + ii * numchan;
      for (jj = 0; jj < numchan; jj++)
          zero_dm_block[ii] += *powptr++;
      ftmp[ii] = zero_dm_block[ii];
   }
   avg_var(ftmp, ptsperblk, &current_avg, &current_std);
   current_std = sqrt(current_std);
   current_med = median(ftmp, ptsperblk);

   /* Calculate the current standard deviation and mean  */
   /* but only for data points that are within a certain */
   /* fraction of the median value.  This removes the    */
   /* really strong RFI from the calculation.            */
   {
      float lo_cutoff, hi_cutoff;
      int numgoodpts = 0;

      lo_cutoff = current_med - 3.0 * current_std;
      hi_cutoff = current_med + 3.0 * current_std;;
      for (jj = 0; jj < numchan; jj++)
         chan_avg_temp[jj] = 0.0;
      /* Find the "good" points */
      for (ii = 0; ii < ptsperblk; ii++) {
         if (zero_dm_block[ii] > lo_cutoff && zero_dm_block[ii] < hi_cutoff) {
            ftmp[numgoodpts] = zero_dm_block[ii];
            powptr = rawdata + ii * numchan;
            for (jj = 0; jj < numchan; jj++)
               chan_avg_temp[jj] += *powptr++;
            numgoodpts++;
         }
      }

      /* Calculate the current average and stddev */
      if (numgoodpts < 1) {
         current_avg = running_avg;
         current_std = running_std;
         for (jj = 0; jj < numchan; jj++)
            chan_avg_temp[jj] = chan_running_avg[jj];
      } else {
         avg_var(ftmp, numgoodpts, &current_avg, &current_std);
         current_std = sqrt(current_std);
         for (jj = 0; jj < numchan; jj++)
            chan_avg_temp[jj] /= numgoodpts;
      }
   }

   /* Update a pseudo running average and stdev */
   if (blocksread) {
      running_avg = 0.9 * running_avg + 0.1 * current_avg;
      running_std = 0.9 * running_std + 0.1 * current_std;
      for (ii = 0; ii < numchan; ii++)
         chan_running_avg[ii] = 0.9 * chan_running_avg[ii] + 0.1 * chan_avg_temp[ii];
   } else {
      running_avg = current_avg;
      running_std = current_std;
      for (ii = 0; ii < numchan; ii++)
         chan_running_avg[ii] = chan_avg_temp[ii];
      if (current_avg == 0.0)
         printf("Warning:  problem with clipping in first block!!!\n\n");
   }

   /* See if any points need clipping */
   trigger = clip_sigma * running_std;
   for (ii = 0; ii < ptsperblk; ii++) {
      if (fabs(zero_dm_block[ii] - running_avg) > trigger) {
         clipit = 1;
         break;
      }
   }

   /* Update the good channel levels */
   for (ii = 0; ii < numchan; ii++)
      good_chan_levels[ii] = chan_running_avg[ii];

   /* Replace the bad channel data with channel median values */
   /* that are scaled to equal the running_avg.               */
   if (clipit) {
      for (ii = 0; ii < ptsperblk; ii++) {
         if (fabs(zero_dm_block[ii] - running_avg) > trigger) {
            powptr = rawdata + ii * numchan;
            for (jj = 0; jj < numchan; jj++)
               *powptr++ = good_chan_levels[jj];
            clipped++;
         }
      }
   }
   blocksread++;

   vect_free(chan_avg_temp);
   vect_free(zero_dm_block);
   vect_free(ftmp);

   return clipped;
}

void apply_mask(float* data, header *h, int64_t nsamples, int64_t offset, float clip_sigma, float *padvals, \
        int *maskchans, mask *obsmask)
// input is the input array of data of the size of nsamples; The type of input is actually dedisp_byte*
// header is metadata, and clip_sigma is the clipping threshold
// padvals is the array of padding values for clipped data
// offset is the offset in samples from the start of the file
// maskchans is an array of length numchans
// contains a list of the number of channels that were masked.
// obsmask is the mask structure to use for masking.
{
   int ii, jj, shift, channum;
   double starttime = offset * h->tsamp, duration = nsamples * h->tsamp;
   int nummasked = 0;
   
   if (obsmask->numchan == 0) return;
   
   /* apply the mask to read input data */
   nummasked = check_mask(starttime, duration, obsmask, maskchans);
           
   /* Clip nasty RFI if requested and we're not masking all the channels */
   if ((clip_sigma > 0.0) && !(nummasked == -1))
      clip_times(data, nsamples, h->nchan, clip_sigma, padvals);
           
   if (nummasked == -1) {     /* If all channels are masked */
       for (ii = 0; ii < nsamples; ii++)
          memcpy(data + ii * h->nchan, padvals, h->nchan * sizeof(float));
   } else if (nummasked > 0) { /* Only some of the channels are masked */
              for (ii = 0; ii < nsamples; ii++) {
                 shift = ii * h->nchan;
                 for (jj = 0; jj < nummasked; jj++) {
                    channum = maskchans[jj];
                    data[shift + channum] = padvals[channum];
                 }
              }
          }
}

/* inverse channel numbers in the mask file to match the order
 * in the filterbank file where first channel is the highest frequency
 * the mask made by rfifind assumed that first channel is low frequency
 */
void inverse_mask(mask *obsmask, float *padvals) {

 for (int ii=0; ii<obsmask->num_zap_chans; ii++) 
     obsmask->zap_chans[ii] = obsmask->numchan - 1 - obsmask->zap_chans[ii];
 for (int ii=0; ii<obsmask->numint; ii++)
     for (int jj=0; jj<obsmask->num_chans_per_int[ii]; jj++)
         obsmask->chans[ii][jj] = obsmask->numchan - 1 - obsmask->chans[ii][jj];
 for (int ii=0; ii<obsmask->numchan; ii++) ELEM_SWAP(padvals[ii], padvals[obsmask->numchan - 1 - ii])
}
