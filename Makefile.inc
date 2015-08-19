PREFIX=/home/bassa/linux
DEDISP=$(PREFIX)/src/dedisp
LOFARDAL=/opt/lofardal
include $(DEDISP)/Makefile.inc

LIBS=-L$(PREFIX)/lib -L$(LOFARDAL)/lib -llofardal -lhdf5 -L$(DEDISP)/lib -ldedisp -lm -lgsl -lgslcblas
INCL=-I$(PREFIX)/include -I. -I$(DEDISP)/include -I$(LOFARDAL)/include
CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -c -g
NVCC_FLAGS=-Xcompiler -Wall -O3
