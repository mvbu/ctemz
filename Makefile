# Automated Makefile
CC = g++
CFLAGS = -g -Wall # -O2
COMPILE = $(CC) $(CFLAGS) -c
OBJFILES := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

LOCAL_LIB_DIR=${HOME}/.local/lib
LOCAL_INC_DIR=${HOME}/.local/include
COMMON_OBJFILES = blzlog.o blzrand.o blzmath.o blzsim.o
CTEMZ_OBJFILES = $(COMMON_OBJFILES) ctemz.o
CTEST_OBJFILES = $(COMMON_OBJFILES) ctest.o
COMMON_HEADER_FILES = blzlog.h blzrand.h blzmath.h blzsim.h

FCOMPILER=gfortran
FCOMPILER_OPTS=-g

all: clean temz ctemz

ctemz: $(CTEMZ_OBJFILES) 
	$(CC) -o ctemz $(CTEMZ_OBJFILES) -L$(LOCAL_LIB_DIR)

%.o: %.cpp
	$(COMPILE) -I$(LOCAL_INC_DIR)/casacore -I$(LOCAL_INC_DIR) -I/usr/local/include/casacore -o $@ $<

clean:
	rm -rf *.o ctemz ctest test

ctest: $(CTEST_OBJFILES) 
	$(CC) -o ctest $(CTEST_OBJFILES) -L$(LOCAL_LIB_DIR)

temz: temz.f
	$(FCOMPILER) $(FCOMPILER_OPTS) -o temz temz.f

temz2: temz2.f
	$(FCOMPILER) $(FCOMPILER_OPTS) -o temz2 temz2.f

# Remove the output files (temz program crashes if you don't do this)
rmo: 
	rm temzspec.txt temzlc.txt temzcheck.txt 

test: test.f
	$(FCOMPILER) $(FCOMPILER_OPTS) -o test test.f

# Copy stuff to Dropbox directory
dropbox:
	cp *.cpp *.f *.h Makefile ~/Dropbox/dev/jet
