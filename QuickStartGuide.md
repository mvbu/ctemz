# Introduction #

ctemz is a C++ port of the Professor Alan Marscher (Boston University) temz code (one file: temz.f) written in Fortran. Both programs model the synchrotron and inverse Compton emissions from blazars (see his [paper](http://arxiv.org/abs/1311.7665)). The temz.f file in this repository is basically his temz.f file from around January 2013 with modifications that I've made to add some purely administrative features such as "fixed" random numbers from a file (described below) and command line arguments.

# Quick Start #

Code should compile and run on a Linux system. After checking out the source code onto your computer:

  * cd to directory containing the code (example: `cd ~/ctemz/trunk`)
  * make a "maps" directory (`mkdir maps`)
  * Compile the `ctemz` executable. The `Makefile` invokes `g++`  (`make ctemz`)
  * Execute the program by typing "`./ctemz -d 1`"  The "1" means to simulate 1 day.

# Command Line Details #

Example command line:

`./ctemz -l info -t -s -d 3 -n 1`

Options:
  * `-l [logLevel]` LOG LEVEL psycho, debug, info, warn, error. Usually set this to warn.
  * `-t [noArgument]` TEST MODE. If this switch is set, run in test mode, which among other things means to use the "fixed" random number mechanism described below
  * -`s [noArgument]` SINGLE THREADED MODE. If this switch is set, run all parallelized code as a single thread, regardless of how many processors a machine has.
  * `-d [daysToSimulate]` SIMULATION TIME. Specify number of days to simulate. Right now this only supports integer days, but it is only a parsing problem that can certainly be fixed.
  * `-n [testOutputID]` TEST OUTPUT. The argument controls which set of test/debug outputs, with current range of 1 through 7 (described in a later section). If this switch is absent or argument is set to 0, no output is generated. Output goes into the file "ctestout.txt" in the current directory.

The output files are written into the current directory and are named: `ctemzlc.txt`, `ctemzspec.txt`, `ctemzpol.txt`, and `maps/ctemzmapxxxx.txt` (where the `xxxx` are decimal digits). Note that `ctemz` will overwrite these files when executed, so save them them to a name or folder if you need to keep the results.

Note that the run-time simulation parameters are read in from a file called `temzinp.txt` which is included in the source distribution. The parameter values can be changed of course, but the file is in source control because it needs to be updated if parameters are added or deleted.

More details on input and output files will be in a different section.