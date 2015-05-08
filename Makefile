HOME_BASSA=/home/bassa/linux
DEDISP=$(HOME_BASSA)/src/dedisp
include $(DEDISP)/Makefile.inc

LIBS=-L$(HOME_BASSA)/lib -llofardal -lhdf5 -L$(DEDISP)/lib -ldedisp -lm
INCL=-I$(HOME_BASSA)/src/DAL -I. -I$(DEDISP)/include
CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -c -g
NVCC_FLAGS=-Xcompiler -Wall -O3
TARGET=dragnet_gpu
all: $(TARGET) strip clean

FILES = $(wildcard *.cu)
OBJS = $(FILES:.cu=.o)

$(TARGET): $(OBJS)
	$(NVCC) $(NVCC_FLAGS) $(OBJS) -o $@ $(LIBS)

%.o : %.cu dragnet_gpu.h
	$(NVCC) $(CFLAGS) $(NVCC_FLAGS) $(INCL) $< -o $@

strip:
	strip $(TARGET)

clean:
	rm -f $(OBJS)
