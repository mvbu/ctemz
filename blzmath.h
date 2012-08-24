#ifndef _INCL_BLZMATH_H_
#define _INCL_BLZMATH_H_

using namespace std;

// Various constants
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
static const double SEC_PER_YEAR = 3.16e7;
static const double SQRT2 = 1.41421;
static const double SQRT3 = 1.73205;
static const double M_PER_PARSEC = 3.086e16;
static const double CM_PER_PARSEC = 3.086e18;
static const double C_CM_PER_SEC = 3e10;
static const double EMC2 = 8.186e-7;
static const double CC2 = 1.29e-9; // TODO: What is this?

// Function pointer type to pass into the Gaussian integration routines
typedef double (*QgFunctionPtr)(double, void *pObject);

// A simple 3-vector type
typedef double BlzVec[3];

class BlzMath {
	// General math routines

 public:
  // Replaces data[0:2*nn-1] by its discrete Fourier transform,  if isign is input as 1
  // or replaces data[0:2*nn-1] by nn times its inverse discrete Fourier transform, if isign is input as -1.
  // data is a double array of length 2*nn, representing a complex array of length nn.
  // nn MUST be an integer power of 2 (this is not checked!) TODO: check for this
	// Based on four1() routine in temz.f
  static void fourier(double data[], int nn, int isign);

  // Returns as ss the integral of the function func between a and  b, 
  // by 5-point Gauss-Legendre integration: the function is evaluated exactly 5
  // times at interior points in the range of integration. Based on qg5() from temz.f
  // which is apparently based on routine from Numerical Recipes
  // pObject is passed to (*pFunction)(). pFunction will know what kind of object it is.
  static double qg5(double a, double b, QgFunctionPtr pFunction, void *pObject);

  // Magnitude from 3 components
  static double mag(double x, double y, double z);
  // Magnitude squared from 3 components
  static double magSquared(double x, double y, double z);
  // 2- and 4-component versions of above
  static double mag(double x, double y);
  static double magSquared(double x, double y);
  static double mag(double x, double y, double z, double w);
  static double magSquared(double x, double y, double z, double w);

  // round the input value (presumably a double or float) to the nearest int
	template<class T>
  static int round(const T input) {
    int low = (int)input;
    int high = (int)(input+1.0);
    T lowDiff = abs(input - (T)low);
    T highDiff = abs(input - (T)high);
    return lowDiff < highDiff ? low : high;
	}

 private:
  BlzMath();
  ~BlzMath();
};
    
#endif // _INCL_BLZMATH_H_

