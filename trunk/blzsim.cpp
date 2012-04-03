#include <cmath>
#include "blzlog.h"
#include "blzmath.h"
#include "blzrand.h"
#include "blzsim.h"

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
