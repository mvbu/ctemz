#include <cmath>
#include "blzlog.h"
#include "blzmath.h"
#include "blzrand.h"
#include "blzsim.h"

BlzSimCommon::BlzSimCommon() {}
BlzSimCommon::~BlzSimCommon() {}
void BlzSimCommon::setGgam(const double gam[CDIST_SIZE]) { for(int i=0; i < CDIST_SIZE; i++)  ggam[i] = gam[i]; }
void BlzSimCommon::setEdist(const double _edist[CDIST_SIZE]) { for(int i=0; i < CDIST_SIZE; i++)  edist[i] = _edist[i]; }

BlzSim::BlzSim() {
}

BlzSim::~BlzSim() {
}

/* Output/result parameter lc_sim must be of dimension BLZSIM_DIM_16384 */
void BlzSim::psdsim(const int N, const float beta1, const float beta2, const float nu_break, const float t_incre1, float *lc_sim)
{
  // This code is ported from Fortran routine temz.f:psdsim() by R. Chatterjee
  double nu[BLZSIM_DIM16384],dat[BLZSIM_DIM16384],R[BLZSIM_DIM32768];
  double dataim[BLZSIM_DIM16384],datareal[BLZSIM_DIM16384],flux_s[BLZSIM_DIM16384];
  double fac_norm,fac_norm2;
  float  ann;
  int j,i,nn,ISEED1,ISEED2,isign;

  fac_norm=1./(N*t_incre1);
  fac_norm2=pow((double)N,2)/(2.*N*t_incre1);

  ISEED1=58;
  ISEED2=256871;
  BlzRand::setIX(ISEED1,ISEED2);
  BlzLog::debugScalarPair("1) ISEED1/ISEED2", ISEED1, ISEED2);
  BlzRand::rnstnr(R,N/2);
  BlzRand::getIX(&ISEED1,&ISEED2);
  BlzLog::debugScalarPair("2) ISEED1/ISEED2", ISEED1, ISEED2);
  for(j=1; j<=N/2; j++) {
    nu[j-1] = j*fac_norm/86400.;
    if(nu[j-1] <= nu_break) {
      // (double) cast necessary on Mac OS X because pow(double, int) not available?? (WTH)
      flux_s[j-1]=::pow(nu[j-1], (double)beta1);
      datareal[j-1]=sqrt((double)(fac_norm2*0.5*flux_s[j-1]))*R[j-1];
    }
    else {
      flux_s[j-1]=(::pow(nu[j-1], (double)beta2)) * ((::pow(nu_break, (double)beta1))/(::pow(nu_break, (double)beta2)));
      datareal[j-1] = sqrt((double)(fac_norm2*0.5*flux_s[j-1]))*R[j-1];
    }
  }
  BlzRand::rnstnr(R,N/2);
  BlzRand::getIX(&ISEED1,&ISEED2);
  BlzLog::debugScalarPair("3) ISEED1/ISEED2", ISEED1, ISEED2);

  for(j=1; j<=N/2; j++) {
    dataim[j-1]=sqrt(0.5*fac_norm2*flux_s[j-1])*R[j-1];
  }

  dat[0]=0.0;
  dat[1]=0.0;
  dat[N]=datareal[(N/2)-1];
  dat[N+1]=0; //  this is at +/-1/T and has to obey A+iB = A-iB => B=0          
  j=1;

  for(i=3; i<=N-1; i+=2) {
    dat[i-1]=datareal[j-1];
    dat[2*(N+1)-i-1]=dat[i-1];
    dat[i]=dataim[j-1];
    dat[2*(N+1)-i] = -1.0*dat[i];
    j=j+1;
  }

  nn=N;
  isign=-1;
  BlzMath::fourier(dat,nn,isign);
  j=1;

  for(i=1; i<=2*nn; i+=2) {
    ann=nn;
    lc_sim[j-1]=dat[i-1]/ann; // taking the real part
    j=j+1;
  }
}

double BlzSim::seedph(const double f)
{
  double sn1, sn2;
  sn1 = common.dsnth1;
  sn2 = common.dsnth2;
  common.freq = f;
  double ans = BlzMath::qg5(sn1, sn2, sdgran, &common);
  // Multiply by 2pi, but then divide by 4pi to get intensity
  double retVal = 0.5 * ans;
  return retVal;
}

// Non-member function
double sdgran(double sn, void *pObject)
{
  BlzSimCommon* pcm = (BlzSimCommon*)pObject;
  const double HOK = 4.80e-11;
  double retVal = 0.0;
  double cs = sqrt((double)(1.0-sn*sn));
  double tdel=1.0/(pcm->gammad*(1.0D-pcm->betad*cs));
  double f = pcm->freq/tdel;
  double expon= HOK * f/pcm->tdust;
  double val = 0.0;

  // The Fortran code has a bunch of goto statements. The structure
  // below implements the same logic, I think :-)
  if(expon <= .01) {
    val=(2.76e-16*f)*pcm->tdust*(f/9.0e20);
  }
  else if(expon > 5.0) {
    if(expon > 15.0)
      return retVal; // seems odd but this is what the Fortran code does
    val=(1.33e-26*f)*f/((9.0e20/f)*exp(expon));
  }
  else {
    val=(1.33e-26*f)*f/((9.0e20/f)*(exp(expon)-1.0));
  }
  retVal = val*tdel*tdel;
  return retVal;
}

double BlzSim::ajnu(const double anu) 
{
  double retVal = 0.0;
  double b = common.bperpp;

  if(b >= 1.0e-5) {
    int i;
    double* gam = common.ggam;
    double* edist = common.edist;
    double a, alg, rfact, gran1, gran2, grat, addit=0.0;
    double xfact = 2.38e-7*anu/b;
    double x = xfact/(gam[0]*gam[0]);
    if(x > .01) {
      rfact = 0.0;
      if(x < 25.0) 
        rfact=(1.08895*::pow(x, (double)0.20949) - 2.35861e-3/::pow(x, (double)0.79051))/exp(x);
    }
    else {
      rfact=1.5*::pow(x, ONETHIRD);
    }

    gran1 = edist[0] * rfact;

    for(i=2; i<=BlzSimCommon::CDIST_SIZE; i++) {
      x = xfact/(gam[i-1]*gam[i-1]);
      if(x > .01) {
        rfact = 0.0;
        if(x < 25.0) 
          rfact=(1.08895*::pow(x, (double)0.20949) - 2.35861e-3/::pow(x, (double)0.79051))/exp(x);
      }
      else {
        rfact=1.5*::pow(x, ONETHIRD);
      }

      gran2 = edist[i-1] * rfact;
      addit = .5 * (gran1+gran2) * (gam[i-1]-gam[i-2]);

      // This block is to implement a bunch of goto logic from the Fortran code
      // Basically, we're deciding if we want the alternate calculation of addit
      if((gran1!=0.0) && (gran2!=0.0)) {
        grat = gam[i-1]/gam[i-2];
        alg = log10(grat);
        if(alg != 0.0) {
          a = 1.0 + log10(gran2/gran1)/alg;
          if((a <= 5.0) && (a >= -5.0)) {
            if(a >= .01)
              addit = gran1 * (::pow(grat,a) - 1.0) * gam[i-2]/a;
          }
        }
      }
    
      retVal = retVal + addit;
      gran1 = gran2;
    } // for(i=2; i<=44; i++) 
  }
  
  return retVal;
}

double BlzSim::akapnu(const double anu)
{
  int i;
  double retVal = 0.0;
  double b = common.bperpp;
  double* gam = common.ggam;
  double* edist = common.edist;
  double xfact = 2.38e-7*anu/b;
  double x = xfact/(gam[0]*gam[0]);
  double f1 = 0.0, f2 = 0.0, addit = 0.0;
  double gran1, gran2, grat, alg, a;

  if(x <= 25.0)
    f1 = 1.8 * ::pow(x, (double).3)/exp(x);

  gran1 = f1 * edist[0]/gam[0];

  for(i=2; i<= BlzSimCommon::CDIST_SIZE; i++) {
    x = xfact / (gam[i-1]*gam[i-1]);
    gran2 = 0.0;
    addit = 0.5 *(gran1+gran2)*(gam[i-1]-gam[i-2]);

    // This block is to implement a bunch of goto logic from the Fortran code
    // Basically, we're deciding if we want the alternate calculation of addit
    if(x <= 25.0) {
      f2 = 1.8 * ::pow(x, (double).3)/exp(x);
      gran2 = f2 * edist[i-1] / gam[i-1];
      // Have to re-calculate addit, since we potentially changed gran2
      addit = 0.5 *(gran1+gran2)*(gam[i-1]-gam[i-2]);

      if((gran1!=0.0) && (gran2!=0.0)) {
        grat = gam[i-1]/gam[i-2];
        alg = log10(grat);
        if(alg != 0.0) {
          a = 1.0 + log10(gran2/gran1)/alg;
          if((a <= 5.0) && (a >= -5.0)) {
            if(a >= .01)
              addit = gran1 * (::pow(grat,a) - 1.0) * gam[i-2]/a;
          }
        }
      }
    } // if(x <= 25.0)

    retVal = retVal + addit;
    gran1 = gran2;
  } // for(i=2; i<= EDIST_SIZE; i++)
  
  return retVal;
}
