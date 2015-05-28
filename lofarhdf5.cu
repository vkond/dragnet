#include "lofarhdf5.h"
#include <iomanip>  // for setprecision use
#include <iostream>

// constructor
HDF5::HDF5(char *filename, header* h, int verbose) : Input() {
  // reading HDF5 input file
  if ((open(filename, h, verbose)) != 0) {
     fprintf(stderr, "Error reading file %s\n", filename);
     exit(-1);
  }
}

// closing the file
void HDF5::close() {
    if (fd != NULL) delete(fd);
    if (stokes != NULL) delete(stokes);
    if (input != NULL) free(input);
}

// destructor
HDF5::~HDF5() {
    close();
}

// opens input HDF5 file and collects all necessary meta data
int HDF5::open(char *filename, header* h, int verbose) {

  current_sample = 0;

  fd = new BF_File (filename);
  cerr << "Reading Lofar HDF5 file..." << endl << endl;
  strcpy(h->inpfile, filename);

  // getting PI
  Attribute<std::string> PI = fd->projectPI();
  if (verbose) if (PI.exists()) cerr << "PI=" << PI.get() << endl;

  // getting project contact
  Attribute<std::string> projectContact = fd->projectContact();
  if (verbose) if (projectContact.exists()) cerr << "PROJECT_CONTACT=" << projectContact.get() << endl;

  // getting targets
  Attribute< std::vector<std::string> > BFtargets = fd->targets();
  if (verbose) {
   if (BFtargets.exists()) {
      std::vector<std::string> t = BFtargets.get();
      std::vector<std::string>::size_type i=0;
      for(std::vector<std::string>::iterator it = t.begin(); it != t.end(); ++it, i++)
        cerr << "TARGET" << i << "=" << *it << endl;
   } else cerr << "TARGET does not exist" << endl;
  }

  // getting frequency center
  Attribute<double> freq = fd->observationFrequencyCenter();
  if (!freq.exists()) {
   cerr << "observationFrequencyCenter not defined" << endl;
   return 1;
  } else if (verbose) cerr << "observation frequency=" << setprecision(20) << freq.get() << endl;

  // getting number of SAPs
  Attribute<unsigned> nsap = fd->observationNofSubArrayPointings();
  if (!nsap.exists()) {
   cerr << "observationNofSubArrayPointings not defined" << endl;
   return 1;
  } else if (verbose) cerr << "number of SAPs=" << nsap.get() << endl;

  // getting  the instance of SAP
  // checking all SAPs if they exist to pick the right one (if there will be two SAPs in one file, only
  // the first one will be picked up)
  unsigned sap_index;
  for (sap_index=0; sap_index<nsap.get(); sap_index++) {
    if (fd->subArrayPointing(sap_index).exists()) break;
  }
  BF_SubArrayPointing sap = fd->subArrayPointing(sap_index);

  Attribute<unsigned> nbeam = sap.observationNofBeams();
  if (!nbeam.exists()) {
   cerr << "sap.observationNofBeams not defined" << endl;
   return 1;
  } else if (verbose) cerr << "number of beams=" << nbeam.get() << endl;

  // getting the instance of first TA beam in the SAP
  // checking all TABs in the SAP if they exist in the file until the first one that exists is found
  unsigned tab_index;
  for (tab_index=0; tab_index<nbeam.get(); tab_index++) {
    if (sap.beam(tab_index).exists()) break;
  }
  BF_BeamGroup beam = sap.beam(tab_index);

  // getting the center frequency of the beam
  Attribute<double> freq2 = beam.beamFrequencyCenter();
  if (!freq2.exists()) {
   cerr << "beam.beamFrequencyCenter not defined" << endl;
   return 1;
  } else if (verbose) cerr << "beam frequency=" << setprecision(20) << freq2.get() << endl;

  // getting the subband width
  Attribute<double> bw2 = beam.subbandWidth();
  if (!bw2.exists()) {
   cerr << "beam.subbandWidth not defined" << endl;
   return 1;
  } else if (verbose) cerr << "sap subbandwidth=" << setprecision(20) << bw2.get() << endl;

  // getting number of channels per subband
  Attribute<unsigned> nchan = beam.channelsPerSubband();
  if (!nchan.exists()) {
   cerr << "beam.channelsPerSubband not defined" << endl;
   return 1;
  } else if (verbose) cerr << "number of channels/sub=" << nchan.get() << endl;

  // getting the pointer for the Stokes class
  stokes = 0;

  for (unsigned i=0; i<4; i++) {
    BF_StokesDataset tmp = beam.stokes(i);
    if (tmp.exists()) {
      stokes = new BF_StokesDataset (beam.stokes(i));
    }
  }

  // getting the Stokes component
  Attribute<std::string> stokesC = stokes->stokesComponent();
  if (verbose) if (stokesC.exists()) cerr << "stokes component=" << stokesC.get() << endl;

  // getting the number of subbands
  Attribute<unsigned> nsub = stokes->nofSubbands();
  if (verbose) { if (nsub.exists()) cerr << "nsub=" << nsub.get() << endl; else cerr << "stokes nofSubbands not defined" << endl; }

  // getting the number of channels for each subband
  Attribute< std::vector<unsigned> > nofchan = stokes->nofChannels();
  if (verbose) if (nchan.exists()) {
                 std::vector<unsigned> nchan = nofchan.get();
                 cerr << "stokes nofChannels size=" << nchan.size() << endl;
                 // for (unsigned i=0; i<nchan.size(); i++) cerr << "stokes nofChannels[" << i << "]=" << nchan[i] << endl;
               } else cerr << "stokes nofChannels not defined" << endl;

  // getting the rank of the dataset
  size_t ndim= stokes->ndims();
  if (verbose) cerr << "stokes ndim=" << ndim << endl;

  if (verbose) {
   std::vector<std::string> files = stokes->externalFiles();
   for (unsigned i=0; i<files.size(); i++)
     cerr << "files[" << i << "]=" << files[i] << endl;
  }

  // getting telescope
  Attribute<std::string> telescope = fd->telescope();
  if (telescope.exists()) {
   if (verbose) cerr << "telescope=" << telescope.get() << endl;
   // setting the telescope
   h->telescope_id = 11;  // For now assuming it is LOFAR 
  }

  // setting machine
  // For now assuming it is LOFAR's COBALT
  h->machine_id = 11;

  // getting the vector of targets
  Attribute< std::vector<std::string> > targets = beam.targets();
  if (targets.exists()) {
    std::vector<std::string> t = targets.get();
    if (t.size() != 0) {
     strcpy(h->source_name, t.front().c_str());
     if (verbose) cerr << "target = " << t.front() << endl;
    } else { if (verbose) cerr << "targets vector is empty" << endl; }
  } else { if (verbose) cerr << "beam target does not exist" << endl; }

  // getting number of samples
  Attribute<unsigned> nsamp = stokes->nofSamples();
  if (nsamp.exists()) h->nsamp = nsamp.get();

  // are data in Complex Voltage format?
  Attribute<bool> volts = beam.complexVoltage();
  if (volts.exists() && volts.get() == 1) {
    cerr << "Can't process complex-voltage data, ndim = " << ndim << endl;
    return 1;
  }
 
  // check for which coordinate is Spectral
  unsigned spectral_dim = 1;

  // getting instance of Coordinates container
  CoordinatesGroup coord = beam.coordinates();
  if (coord.exists()) {
    Attribute< std::vector<std::string> > types = coord.coordinateTypes();
    if (types.exists()) {
      std::vector<std::string> t = types.get();
      for (unsigned i=0; i<t.size(); i++) {
	if (t[i] == "Spectral") {
	  spectral_dim = i;
	  break;
	}
      }
    }
  }

  std::vector<ssize_t> dims = stokes->dims();
  h->nchan = dims[spectral_dim];
  cerr << "Total number of channels=" << h->nchan << endl;
  
  // getting number of Stokes components in one file
  //Attribute<unsigned> npol = beam.nofStokes();
  // getting number of Stokes components in the observation
  Attribute<unsigned> npol = beam.observationNofStokes();
  unsigned stokes_npol = 1;

  if (npol.exists()) stokes_npol = npol.get();

  if (stokes_npol == 1) {
    h->nif = 1;
  } else {
    cerr << "Can't process more than one IFs" << endl;
    return 1;
  }

  h->nbit = 32;
  h->nbeam = 1; // For now assuming there is only one beam
  h->ibeam = 0;

  // getting split Frequency center of the beam
  Attribute<double> cfreq = beam.beamFrequencyCenter();
  if (!cfreq.exists()) {
   cerr << "beamFrequencyCenter not defined" << endl;
   return 1;
  } else { if (verbose) cerr << "beamFrequencyCenter=" << setprecision(20) << cfreq.get() << endl;
  }

  // getting the start MJD
  Attribute<double> mjd = fd->observationStartMJD();
  if (mjd.exists()) h->tstart = mjd.get();
  if (verbose) cerr << "MJD=" << setprecision(20) << h->tstart << endl;

  // getting the clock rate
  Attribute<double> cRate = fd->clockFrequency();
  if (verbose) {
    if (cRate.exists()) cerr << "clockRate=" << setprecision(20) << cRate.get() << endl;
    else cerr << "clockRate undefined" << endl;
  }

  // getting the sampling rate
  Attribute<double> sRate = beam.samplingRate();
  if (verbose) {
    if (sRate.exists()) cerr << "samplingRate=" << setprecision(20) << sRate.get() << endl;
    else cerr << "samplingRate undefined" << endl;
  }

  // getting the sampling time
  Attribute<double> sTime = beam.samplingTime();
  if (sTime.exists()) {
    if (verbose) cerr << "samplingTime=" << setprecision(20) << sTime.get() << " s"<< endl;
    h->tsamp = sTime.get();
  } else if (verbose) cerr << "samplingTime undefined" << endl;

  // getting the channel width
  Attribute<double> rate = beam.channelWidth();
  if (!rate.exists()) {
   cerr << "beam.channelWidth not defined" << endl;
   return 1;
  } else { if (verbose) cerr << "channel Width=" << setprecision(20) << rate.get() << " Hz" << endl;
           h->foff = -1. * rate.get() * 1.0e-6; // We make it negative as is needed by Sigproc, and in MHz !!
         }

  // getting the subband width
  Attribute<double> subwidth = beam.subbandWidth();
  if (verbose) if (subwidth.exists())
                 cerr << "subband Width=" << setprecision(20) << subwidth.get() << " Hz" << endl;
               else cerr << "subband Width undefined" << endl;

  // setting the bandwidth (in MHz) of the file
  double bw_file = h->nchan * rate.get() * 1.0e-6;
  h->fch1 = cfreq.get() + bw_file/2. - fabs(h->foff * 1.0e-6)/2.; // in MHz !

  // getting the RA and DEC of the beam (in degrees)
  Attribute<double> radeg = beam.pointRA();
  if (verbose) {
    if (radeg.exists()) {
      cerr << "RA=" << setprecision(20) << radeg.get() << " deg" << endl;
      int ra_h = (int)(radeg.get()/15.);
      int ra_m = (int)((radeg.get()/15. - ra_h)*60.);
      double ra_s = (radeg.get()/15. - ra_h - ra_m/60.)*3600.;
      char tmp[30];
      sprintf(tmp, "%02d%02d%s%lf", ra_h, ra_m, ra_s < 10 ? "0" : "", ra_s);
      sscanf(tmp, "%lf", &h->src_raj);
      cerr << "RA=" << tmp << endl;
    } else cerr << "RA undefined" << endl;
  }

  Attribute<double> decdeg = beam.pointDEC();
  if (verbose) {
    if (decdeg.exists()) {
      cerr << "DEC=" << setprecision(20) << decdeg.get() << " deg" << endl;
      int dec_d = (int)(fabs(decdeg.get()));
      int dec_m = (int)((fabs(decdeg.get()) - dec_d)*60.);
      double dec_s = (fabs(decdeg.get()) - dec_d - dec_m/60.)*3600.;
      int sign = (int)(decdeg.get());
      if (sign < 0) dec_d = -dec_d;
      char tmp[30];
      sprintf(tmp, "%02d%02d%s%lf", dec_d, dec_m, dec_s < 10 ? "0" : "", dec_s);
      sscanf(tmp, "%lf", &h->src_dej);
      cerr << "DEC=" << tmp << endl;
    } else cerr << "DEC undefined" << endl;
  }

  return 0;
}


// reading the nsamples from the file
int64_t HDF5::read(int64_t nsamples, int64_t shift_back, header* h, void*& out) {

  // how many samples we actually read
  int64_t read_samples = current_sample + nsamples > h->nsamp ? h->nsamp - current_sample : nsamples;
  // allocating memory for input buffer
  input=(dedisp_float *) realloc(input, sizeof(dedisp_float) * read_samples * h->nchan);
  // Reading the data
  vector<size_t> pos (2);
  pos[0] = current_sample;
  pos[1] = 0;

  // I need this because freq-order is different: lowest freq is first in LOFAR data
  // and for dedisp library it's required the highest freq to be the first
  // so we need to re-order
  float *outbuf = new float[read_samples * h->nchan];
  stokes->get2D (pos, outbuf, read_samples, h->nchan);
  for (long ii = 0; ii < read_samples; ii++) {
   for (long jj = 0; jj < h->nchan; jj++) {
    input[(ii+1)*h->nchan-1-jj] = outbuf[ii*h->nchan+jj];
//    cerr << "samp=" << ii << "  chan=" << jj << "  val=" << outbuf[(ii+1)*h->nchan-1-jj] << endl;
   }
  }
  delete(outbuf);

  out = input; // set the (void*) pointer to where the data are
  current_sample += (read_samples - shift_back);
  return read_samples - shift_back;
}
