# Automated Makefile
CC = g++
CFLAGS = -g -Wall # -O2
COMPILE = $(CC) $(CFLAGS) -c
OBJFILES := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

LOCAL_LIB_DIR=${HOME}/.local/lib
LOCAL_INC_DIR=${HOME}/.local/include
COMMON_OBJFILES = blzlog.o blzrand.o blzmath.o blzsim.o blzsiminputreader.o
CTEMZ_OBJFILES = $(COMMON_OBJFILES) ctemz.o
CTEST_OBJFILES = $(COMMON_OBJFILES) ctest.o
SCRATCH_OBJFILES = $(COMMON_OBJFILES) scratch.o
COMMON_HEADER_FILES = blzlog.h blzrand.h blzmath.h blzsim.h blzsiminput.h blzsiminputreader.h

FCOMPILER=gfortran
FCOMPILER_OPTS=-g -ffixed-line-length-none

all: clean temz ctemz scratch

ctemz: $(CTEMZ_OBJFILES) 
	$(CC) -o ctemz $(CTEMZ_OBJFILES) -L$(LOCAL_LIB_DIR)

%.o: %.cpp
	$(COMPILE) -I$(LOCAL_INC_DIR) -o $@ $<

clean:
	rm -rf *.o ctemz ctest test scratch

ctest: $(CTEST_OBJFILES) 
	$(CC) -o ctest $(CTEST_OBJFILES) -L$(LOCAL_LIB_DIR)

scratch: $(SCRATCH_OBJFILES) 
	$(CC) -o scratch $(SCRATCH_OBJFILES) -L$(LOCAL_LIB_DIR)

temz: temz.f
	$(FCOMPILER) $(FCOMPILER_OPTS) -o temz temz.f

# Remove the output files (temz.f program crashes if you don't do this)
rmo: 
	rm temzspec.txt temzlc.txt temzpol.txt

test: test.f
	$(FCOMPILER) $(FCOMPILER_OPTS) -o test test.f

randtest: randtest.f
	$(FCOMPILER) $(FCOMPILER_OPTS) -o randtest randtest.f

# Copy stuff to Dropbox directory
dropbox:
	cp *.cpp *.f *.h Makefile ~/Dropbox/dev/jet
