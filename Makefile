# Automated Makefile
CC = g++ -fopenmp
CFLAGS = -g -mcmodel=medium # -Wall -O2
COMPILE = $(CC) $(CFLAGS) -c
OBJFILES := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

LOCAL_LIB_DIR=${HOME}/.local/lib
LOCAL_INC_DIR=${HOME}/.local/include
MPI_LIB_DIR=/usr/lib64/openmpi/lib
MPI_INC_DIR=/usr/include/openmpi-x86_64
COMMON_OBJFILES = blzutil.o blzlog.o blzrand.o blzmath.o blzsiminputreader.o
CTEMZ_OBJFILES = $(COMMON_OBJFILES) ctemz.o
CTEST_OBJFILES = $(COMMON_OBJFILES) ctest.o
CECDUST_OBJFILES = $(COMMON_OBJFILES) cecdust.o
CSSC_OBJFILES = $(COMMON_OBJFILES) cssc.o
SCRATCH_OBJFILES = scratch.o blzutil.o
MPITEST_OBJFILES = mpitest.o
GOMPTEST_OBJFILES = gomptest.o
SERTEST_OBJFILES = sertest.o
COMMON_HEADER_FILES = blzlog.h blzrand.h blzmath.h blzsim.h blzsiminput.h blzsiminputreader.h

FCOMPILER=gfortran
FCOMPILER_OPTS=-g -ffixed-line-length-none

all: clean temz ctemz scratch

ctemz: $(CTEMZ_OBJFILES) blzsim.o
	$(CC) -o ctemz $(CTEMZ_OBJFILES) blzsim.o -L$(LOCAL_LIB_DIR)

ctemzp: $(CTEMZ_OBJFILES) blzsimp.o
	$(CC) -o ctemzp $(CTEMZ_OBJFILES) blzsimp.o -L$(LOCAL_LIB_DIR)

%.o: %.cpp
	$(COMPILE) -I$(LOCAL_INC_DIR) -I$(MPI_INC_DIR) -o $@ $<

clean:
	rm -rf *.o ctemz ctest test scratch

ctest: $(CTEST_OBJFILES) 
	$(CC) -o ctest $(CTEST_OBJFILES) -L$(LOCAL_LIB_DIR)

cecdust: $(CECDUST_OBJFILES) 
	$(CC) -o cecdust $(CECDUST_OBJFILES) -L$(LOCAL_LIB_DIR)

cssc: $(CSSC_OBJFILES) 
	$(CC) -o cssc $(CSSC_OBJFILES) -L$(LOCAL_LIB_DIR)

scratch: $(SCRATCH_OBJFILES) 
	$(CC) -o scratch $(SCRATCH_OBJFILES) -L$(LOCAL_LIB_DIR)

mpitest: $(MPITEST_OBJFILES) 
	$(CC) -o mpitest $(MPITEST_OBJFILES) -L$(MPI_LIB_DIR) -L$(LOCAL_LIB_DIR) -lmpi -lmpi_cxx

gomptest: $(GOMPTEST_OBJFILES) 
	$(CC) -o gomptest $(GOMPTEST_OBJFILES) -L$(LOCAL_LIB_DIR)

sertest: $(SERTEST_OBJFILES) 
	$(CC) -o sertest $(SERTEST_OBJFILES) -L$(LOCAL_LIB_DIR)

temz: temz.f
	$(FCOMPILER) $(FCOMPILER_OPTS) -o temz temz.f

temzd: temzd8f.f
	$(FCOMPILER) $(FCOMPILER_OPTS) -o temzd temzd8f.f

# Remove the output files (temz.f program crashes if you don't do this)
rmo: 
	rm temzspec.txt temzlc.txt temzpol.txt
rmod: 
	rm temzdspec.txt temzdlc.txt temzdpol.txt temzdcheck.txt

test: test.f
	$(FCOMPILER) $(FCOMPILER_OPTS) -o test test.f

ecdust: ecdust.f
	$(FCOMPILER) $(FCOMPILER_OPTS) -o ecdust ecdust.f

ssc: ssc.f
	$(FCOMPILER) $(FCOMPILER_OPTS) -o ssc ssc.f

randtest: randtest.f
	$(FCOMPILER) $(FCOMPILER_OPTS) -o randtest randtest.f

# Work in progress to Dropbox directory. Should then be synced to the remote Dropbox by Dropbox daemon.
dropbox:
	cp *.cpp *.f *.h temz*.txt ctemz*.txt Makefile ~/Dropbox/dev/jet

# Copy the CVS repository to Dropbox directory. Should then be synced to the remote Dropbox by Dropbox daemon.
cvsbackup:
	cp -f -r ~/CVSROOT ~/Dropbox/dev
