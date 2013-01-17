#include <string>
#include <sstream>
#include <cstring>
#include <cmath>
#include "blzmath.h"

void BlzMath::fourier(double data[], int nn, int isign)
{
  //  Replaces data(1:2*nn) by its discrete Fourier transform, 
  //  if isign is input as 1; or replaces
  //  data(1:2*nn) by nn times its inverse discrete Fourier 
  //  transform, if isign is input as -1.
  //  data is a complex array of length nn or, equivalently, 
  //  a real array of length 2*nn.
  //  nn MUST be an integer power of 2 (this is not checked!).
  int i,istep,j,m,mmax,n;
  double tempi,tempr;
  double theta,wi,wpi,wpr,wr,wtemp; // Double precision for the trigonometric recurrences. 
  n=2*nn;
  j=1;

  for(i=1; i<=n; i+=2) {
    if(j > i) {
      tempr=data[j-1]; // Exchange the two complex numbers.
      tempi=data[j];
      data[j-1]=data[i-1];
      data[j]=data[i];
      data[i-1]=tempr;
      data[i]=tempi;
    }
    m=nn;

    while((m>=2) && (j>m)) {
      j=j-m;
      m=m/2;
    }

    j=j+m;
  }

  // Danielson-Lanczos section of the routine
  mmax=2;
  // outer loop executed log2 nn times
  while(n>mmax) {
    istep = 2 * mmax;
    theta=TWOPI/(isign*mmax); //Initialize for the trigonometric recurrence
    wpr=double(-2)*pow(sin((double)0.5*theta),2);
    wpi=sin(theta);
    wr=1.e0;
    wi=0.e0;

    for(m=1; m<=mmax; m+=2) {
      for(i=m; i<=n; i+=istep) {
        j=i+mmax;  //This is the Danielson-Lanczos formula:
        tempr=float(wr)*data[j-1]-float(wi)*data[j];
        tempi=float(wr)*data[j]+float(wi)*data[j-1];
        data[j-1]=data[i-1]-tempr;
        data[j]=data[i]-tempi;
        data[i-1]=data[i-1]+tempr;
        data[i]=data[i]+tempi;
      }
      wtemp=wr; // trig recurrence
      wr=wr*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }

    mmax=istep;
  } // end while(n > mmax)
}

double BlzMath::qg10(double a, double b, QgFunctionPtr pFunction, void *pObject) 
{
  int j;
  // The abscissas and weights:
  double dx,xm,xr;
  double w[5] = {.2955242247, .2692667193, .2190863625, .1494513491, .0666713443};
  double x[5] = {.1488743389, .4333953941, .6794095682, .8650633666, .9739065285};
  xm = 0.5*(b+a);
  xr = 0.5*(b-a);
  double ss = 0.0;
  // Will be twice the average value of the function, since the ten
  // weights (five numbers above each used twice) sum to 2. 
  for(j=1; j<=5; j++) { // do 11 j=1,5
    double dx = xr*x[j-1];
    ss = ss+w[j-1]*((*pFunction)(xm+dx, pObject)+(*pFunction)(xm-dx, pObject));
  } //11 continue

  // Scale the answer to the range of integration.
  ss = xr*ss;
  return ss;
}

// Returns the integral of the function func between a and  b, by 5-point Gauss-Legendre integration:
// the function is evaluated exactly five times at interior points in the range of integration.
// pObject is passed to (*pFunction)()
double BlzMath::qg5(double a, double b, QgFunctionPtr pFunction, void *pObject) 
{
  double retVal = 0.0;
  int j;
  double xm, xr;
  // save w,x <--The Fortran version does this, meaning values in w and x should be retained
  // between calls. But w and x values never get changed from the initial values below, so until
  // I figure out a/the reason why they have to be saved, I'm leaving that part out. -msv April2012
  // The abscissas and weights:
  double w[5] = {.236926885,.478628670,.568888889,.47862867, .236926885 };
  double x[5] = {-.906179846,-.538469310,0.0,.53846931,.906179846 };
  xm=0.5*(a+b);
  xr=0.5*(b-a);

  for(j=0; j<5; j++) {
    retVal=retVal+w[j]*(*pFunction)(xm+x[j]*xr, pObject);
  }

  // Scale the answer to the range of integration.
  retVal=xr*retVal;
  return retVal;
}

// Rotate unit vector a by an angle psi about unit vector c along a great circle to create new vector v
void BlzMath::vecRot(const double ax, const double ay, const double az, // a: unit vector to rotate
                     const double cx, const double cy, const double cz, // c: unit vector rotation axis
                     const double psi,                                  // rotation angle about c
                     double* vx, double* vy, double* vz)                // v: resulting vector
{
  double cs = cos(psi);
  double s = sin(psi);
  double cs1 = 1.0 - cs;
  *vx = ax*(1.0+cs1*(cx*cx-1.0))-ay*(cz*s-cs1*cx*cy)+az*(cy*s+cs1*cx*cz);
  *vy = ax*(cz*s+cs1*cx*cy)+ay*(1.0+cs1*(cy*cy-1.0))-az*(cx*s-cs1*cy*cz);
  *vz = ax*(-cy*s+cs1*cx*cz)+ay*(cx*s+cs1*cy*cz)+az*(1.0+cs1*(cz*cz-1.0));
}

double BlzMath::mag(double x, double y, double z) {  return sqrt(magSquared(x,y,z)); }
double BlzMath::magSquared(double x, double y, double z){ return x*x + y*y + z*z; }
double BlzMath::mag(double x, double y) {  return sqrt(magSquared(x,y)); }
double BlzMath::magSquared(double x, double y){ return x*x + y*y; }
double BlzMath::mag(double x, double y, double z, double w) {  return sqrt(magSquared(x,y,z, w)); }
double BlzMath::magSquared(double x, double y, double z, double w){ return x*x + y*y + z*z + w*w; }



