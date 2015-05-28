#include "dragnet.h"
#include "input.h"

using namespace std;
#include "dal/lofar/BF_File.h"
using namespace dal;


// class to read input data in HDF5 format 
class HDF5 : public Input {

    BF_File* fd; // file descriptor
    BF_StokesDataset* stokes;
    dedisp_float* input;

    public:
        HDF5(char *filename, header* h, int verbose=1);
        ~HDF5(); 

        int open(char *filename, header* h, int verbose=1);
        int64_t read(int64_t nsamples, int64_t shift_back, header* h, void*& out); // read number of samples
        void close();
};
