#ifndef _INCL_BLZSIM_H_
#define _INCL_BLZSIM_H_

using namespace std;

// I'll give these better names when I know what they mean
static const int BLZSIM_DIM16384 = 16384;
static const int BLZSIM_DIM32768 = 32768;

/**
   Various simulation routines
 */
class BlzSim {

 public:
  // code is ported from Fortran psdsim() in temz.f from Ritaban Chatterjee
  // to create variations according to an input PSD with slopes beta1 and beta2 at
  // variational frequencies below and above a break frequency nu_break */
  // N=number of data points in the lc, should be an integer power of 2, N<=8192
  // t_incre1=increment in time in each step while resampling the simulated data.
  // This must be larger than the smallest interval between successive data points in the input light curve
  static void psdsim(const int N, const float beta1, const float beta2, const float nu_break, const float t_incre1, 
                     float *lc_sim);

 private:
  BlzSim();
  ~BlzSim();
};
    
#endif // _INCL_BLZSIM_H_

