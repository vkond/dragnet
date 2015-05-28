#include "dragnet.h"

// class to read filterbank input data in SIGPROC format 
class Sigproc : public Input {

    FILE* fd; // file descriptor
    dedisp_byte* input;

    int read_sigproc_header(FILE* file, header* h); // reading Sigproc header

    public:
        Sigproc(char *filename, header* h, int verbose=1);
        ~Sigproc(); 
    protected:
        int open(char *filename, header* h, int verbose=1);
        int64_t read(int64_t nsamples, int64_t shift_back, header* h, void*& out); // read number of samples
        void close();
};
