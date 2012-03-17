#ifndef _INCL_BLZMATH_H_
#define _INCL_BLZMATH_H_

using namespace std;

static const double TWOPI =     6.28318530717959;
static const double PI =        3.14159265358979;
static const double PIOVERTWO = 1.57079632679490;
static const double RAD_PER_DEG = PI/180.;
static const double DEG_PER_RAD = 180./PI;
static const double SEC_PER_HOUR = 3600.;
static const double SEC_PER_DAY = 86400.;
static const double MIN_PER_HOUR = 60.;
static const double MIN_PER_DAY = 1440.;
static const double HOUR_PER_DAY = 24.;

class BlzMath {
	// General math routines

 public:
	// Based on four1() routine in temz.f
  static void fourier(double data[], int nn, int isign);

 private:
  BlzMath();
  ~BlzMath();
};
    
#endif // _INCL_BLZMATH_H_

