# Fortran Version #

As of the genesis of this Google Code project (February 2014), the `ctemz` code in this repository is a C++ port of the `temz.f` code as it existed in January 2013. Some purely administrative changes were made to `temz.f` as described below.

# C++ Porting notes #

There are certainly some "gotchas" when porting from Fortran to C++, and vice-versa. Here are a few of them in case any of it is helpful:

  * C++ arrays are zero-based, while Fortran arrays are one-based
  * When casting to an integer, C++ truncates while Fortran rounds up or down.



# Changes made to `temz.f` #

Changes were made to `temz.f` to facilitate the C++ port itself and testing of the C++ code, namely, comparing its output with the output of the Fortran code (`temz.f`).  It should be emphasized that no changes were made to any of the physics/modeling part of the code - changes were purely administrative.

## Fixed Random Numbers ##
The most significant change was to implement the Fixed Random Number (FRN) generator as had been implemented in the C++ port, essentially doing some porting in the "reverse" direction. The FRN is a system that reads "random" numbers from a text file containing pre-generated, presumably random, numbers. FRN was essential for the porting process so that both `temz` and `ctemz` would be using the exact same numbers and so the output of both executables could be compared directly to see if the port (the C++ code) produced the same output numbers that the original (Fortran code) produced. Both executables will use the FRN system if the `-t` command line switch is used. The file that contains the carriage-return-separated numbers must be named `fixedrand.dat` and contain 100,000 numbers. Certainly both of these requirements can be made less restrictive in a future release.


## Command Line Argument support ##
The Fortran code was modified to support the same command line arguments as the C++ code. The one exception to this is the `-s` option which tells the C++ code to not do any multi-threading (that is, to not run anything in parallel). No parallelization has been implemented in the Fortran code, so the `-s` option is not supported in the Fortran code (i.e., in `temz`).

## Test Output Modes ##
The same test output modes present in the C++ code have been added to the Fortran code to make it easier to verify the C++ port and/or debug the C++ port when the outputs do not match. The various test output modes are described in a separate section.