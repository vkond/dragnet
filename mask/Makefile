include ../Makefile.inc

CFLAGS_EXTRA = -fpermissive

HFILES = $(wildcard *.h)
CFILES = $(wildcard *.c)
OBJS = $(CFILES:.c=.o)

# target name here should be the same as directory name
mask: $(OBJS)
	@echo

%.o : %.c $(HFILES)
	$(GXX) $(CFLAGS) $(CFLAGS_EXTRA) $(INCL) $< -o $@

clean:
	rm -f $(OBJS)
