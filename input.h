#include "dragnet.h"

#ifndef _INPUT
#define _INPUT

// interface class to open/read input data
class Input {
    public:
        int64_t current_sample;

        Input() {}; // constructor
        ~Input() {}; // destructor

        virtual int open(char *filename, header* h, int verbose=1) = 0;
        virtual int64_t read(int64_t nsamples, int64_t shift_back, header* h, void*& out) = 0; // read number of samples
        virtual void close() = 0;
};

#endif
