PREFIX=/home/vlad/pulsar
DEDISP=$(PREFIX)/src/dedisp-multi
include $(DEDISP)/Makefile.inc

LIBS=-L$(PREFIX)/lib -llofardal -lhdf5 -L$(DEDISP)/lib -ldedisp -lm
INCL=-I$(PREFIX)/include -I. -I$(DEDISP)/include
CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -c -g
NVCC_FLAGS=-Xcompiler -Wall -O3
