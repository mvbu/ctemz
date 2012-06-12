#ifndef _INCL_BLZSIM_H_
#define _INCL_BLZSIM_H_

#include "blzsiminput.h"
#include "blzrand.h"

using namespace std;

// I'll give these better names when I know what they mean
static const int BLZSIM_DIM16384 = 16384;
static const int BLZSIM_DIM32768 = 32768;

// These are variables in the TEMZ (Fortran) "common" blocks. Putting
// them here in a structure for now, to expedite porting to C++.
// Can come up with a better way later.
class BlzSimCommon {
 public:
  static const int CDIST_SIZE = 44;
  static const int CSEED_SIZE = 22;
  static const int CSSC_SIZE = 68;

  double bdx, bdy, bdz, gamd, betd; // cvel
  double zred1, bfld, bperpp; // cparm
  double dcsth1, dcsth2, dsnth1, dsnth2, dsang, tdust; // cdust
  // The first two entries in the next two lines are both referred to in the Fortran as dnu, di
  double snu[CSSC_SIZE], ssseed[CSSC_SIZE]; // cssc 
  double dustnu[CSEED_SIZE],dusti[CSEED_SIZE];    // cseed
  int    nuhi; // cssc
  double ggam[CDIST_SIZE], edist[CDIST_SIZE]; // cdist
  double freq; // cfreq

  BlzSimCommon();
  ~BlzSimCommon();
  // Setters for the arrays. For the non-array members, just set directly for now
  void setGgam(const double gam[CDIST_SIZE]);
  void setEdist(const double _edist[CDIST_SIZE]);
  void setDustnu(const double dnu[CSEED_SIZE]);
  void setDusti(const double di[CSEED_SIZE]);
  void setSnu(const double _snu[CSSC_SIZE]);
  void setSsseed(const double _ssseed[CSSC_SIZE]);
};


// Callback function passed to BlzMath::qg5()
double sdgran(double sn, void *pObject);

/**
   Top-level simulation object and supporting methods
 */
class BlzSim {

 public:
  BlzSim();
  ~BlzSim();
  
  // This is the main method of this class. Eventually, all other methods will be private.
  // Instantiate a BlzSim object, instantiate and initialize a BlzSimInput object, then pass it
  // to the run() method. ndays is number of days to simulate
  void run(BlzSimInput& blzSimInput, double ndays=20.0, bool bTestMode = false);

  // Most of these methods will eventually become private or protected.

  // Instance of structure to store shared variables, serving same purpose
  // as the COMMON variables in the Fortran code
  BlzSimCommon common;
  // Instance of the BlzRand object common to this sim instance
  BlzRand randObj;

  // code is ported from Fortran psdsim() in temz.f from Ritaban Chatterjee
  // to create variations according to an input PSD with slopes beta1 and beta2 at
  // variational frequencies below and above a break frequency nu_break
  // N=number of data points in the lc, should be an integer power of 2, N<=8192
  // t_incre1=increment in time in each step while resampling the simulated data.
  // This must be larger than the smallest interval between successive data points in the input light curve
  void psdsim(const int N, const double beta1, const double beta2, const double nu_break, 
              const double t_incre1, double *lc_sim);

  double seedph(const double f);
  double ajnu(const double anu);
  double akapnu(const double anu);
  double ecdust(const double anuf);
  double ssc(const double anuf);
  double polcalc(const double b, const double bx, const double by, const double bz, 
                 const double clos, const double slos);
 private:
  void initRandFromTime(bool bTestMode = false);
  static const double ONETHIRD = .33333333;
  static const double S0 = 6.237;
};
    
#endif // _INCL_BLZSIM_H_

