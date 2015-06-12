include Makefile.inc

# name of the output binary (i.e. project)
TARGET=dragnet
# sub-directories (should be separated by space)
DIRS = mask
all: $(TARGET) strip clean

DIROBJS = $(DIRS:%=%/*.o)

HFILES = $(wildcard *.h)
CFILES = $(wildcard *.cxx)
OBJS = $(CFILES:.cxx=.o)

CUFILES = $(wildcard *.cu)
CUDAOBJS = $(CUFILES:.cu=.obj)

$(TARGET): $(OBJS) $(CUDAOBJS) $(DIRS)
	$(NVCC) $(NVCC_FLAGS) $(OBJS) $(CUDAOBJS) $(DIROBJS) -o $@ $(LIBS)

%.o : %.cxx $(HFILES)
	$(GXX) $(CFLAGS) $(INCL) $< -o $@

%.obj : %.cu $(HFILES)
	$(NVCC) $(NVCC_FLAGS) $(CFLAGS) $(INCL) $< -o $@

$(DIRS) : .PHONY
	@cd $@; $(MAKE) $@ ; cd ..

strip:
	strip $(TARGET)

clean:
	rm -f $(OBJS) $(CUDAOBJS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

.PHONY :
