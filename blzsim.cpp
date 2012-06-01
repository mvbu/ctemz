#include <cmath>
#include "blzlog.h"
#include "blzmath.h"
#include "blzrand.h"
#include "blzsim.h"

BlzSimCommon::BlzSimCommon() {}
BlzSimCommon::~BlzSimCommon() {}
void BlzSimCommon::setGgam(const double gam[CDIST_SIZE]) { for(int i=0; i < CDIST_SIZE; i++)  ggam[i] = gam[i]; }
void BlzSimCommon::setEdist(const double _edist[CDIST_SIZE]) { for(int i=0; i < CDIST_SIZE; i++)  edist[i] = _edist[i]; }
void BlzSimCommon::setDustnu(const double dnu[CSEED_SIZE]) { for(int i=0; i < CSEED_SIZE; i++)  dustnu[i] = dnu[i]; }
void BlzSimCommon::setDusti(const double di[CSEED_SIZE]) { for(int i=0; i < CSEED_SIZE; i++)  dusti[i] = di[i]; }
void BlzSimCommon::setSnu(const double _snu[CSSC_SIZE]) { for(int i=0; i < CSSC_SIZE; i++)  snu[i] = _snu[i]; }
void BlzSimCommon::setSsseed(const double _ssseed[CSSC_SIZE]) { for(int i=0; i < CSSC_SIZE; i++)  ssseed[i] = _ssseed[i]; }

BlzSim::BlzSim() {
}

BlzSim::~BlzSim() {
}

/* Output/result parameter lc_sim must be of dimension BLZSIM_DIM_16384 */
void BlzSim::psdsim(const int N, const float beta1, const float beta2, const float nu_break, const float t_incre1, float *lc_sim)
{
  // This code is ported from Fortran routine temz.f:psdsim() by R. Chatterjee
  double nu[BLZSIM_DIM16384],dat[BLZSIM_DIM32768],R[BLZSIM_DIM32768];
  double dataim[BLZSIM_DIM16384],datareal[BLZSIM_DIM16384],flux_s[BLZSIM_DIM16384];
  double fac_norm,fac_norm2;
  float  ann;
  int j,i,nn,ISEED1,ISEED2,isign;

  fac_norm=1./(N*t_incre1);
  fac_norm2=pow((double)N,2)/(2.*N*t_incre1);

  // randObj is an instance of BlzRand inside this BlzSim instance
  ISEED1=58;
  ISEED2=256871;
  randObj.setIX(ISEED1,ISEED2);
  BlzLog::debugScalarPair("1) ISEED1/ISEED2", ISEED1, ISEED2);
  randObj.rnstnr(R,N/2);
  randObj.getIX(&ISEED1,&ISEED2);
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
  randObj.rnstnr(R,N/2);
  randObj.getIX(&ISEED1,&ISEED2);
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

double BlzSim::ecdust(const double anuf)
{
  // This function (like its original Fortran counterpart) contains two loops that 
  // are nearly identical. Hence this C++ version has the same "flaw" just to port it more
  // quickly. Eventually can write one function that is called for both loops. Hopefully 
  // this function is general enough to be called for both loops in ecdust(), as well.
  double* gam = common.ggam;
  double* edist = common.edist;
  double* dnu = common.dustnu;
  double* di = common.dusti;
  double retVal = 0.0, gran = 0.0;
  // (MSV: this is from the Fortran and the line of code is 
  // commented out in the Fortran, so including it all here for reference)
  // In this version, approximate that plasma velocity is along jet axis
  // and that Doppler factor of dust torus is the mean over its solid angle
  // as viewed in the plasma frame 
  // tdel=1.0/(gammad*(1.0d0-betad*dcsth1))
  double g1=gam[0];
  double vala=S0*edist[0]/(g1*g1);
  int ie;

  //
  // Top-level loop to integrate over electron Lorentz factors (gam)
  //
  for(ie=1; ie<=BlzSimCommon::CDIST_SIZE-1; ie++) {
    double g2 = gam[ie], die;
    double valb = S0 * edist[ie]/(g2*g2);
    double val1=0.0, val2=0.0, gran1=0.0, addit=0.0;
    // Set up loop 1 to integrate over incident photon frequency anui
    double anumax=min(anuf,dnu[BlzSimCommon::CSEED_SIZE-1]);
    double di1=di[0];
    int id=2;
    double anumin=0.25*anuf/(g1*g1);
    // skip over loop 1 entirely if this condition is not satisfied
    if(anumin <= anumax) {
      if(anumin <= dnu[0]) {
        anumin = dnu[0];
      }
      else {
        for(id=2; id<=BlzSimCommon::CSEED_SIZE; id++) {
          if(anumin <= dnu[id-1])
            break;
        }
        // if (anumin <= dnu[id-1]) is satisfied before id=CSEED_SIZE, then id will be set
        // to that number. Otherwise it will get to CSEED_SIZE (currently 22). The Fortran
        // code explicitly sets id to 22 in this case, but I don't see why - id will already
        // have the value 22
        double a = log10(di[id-2]/di[id-1])/log10(dnu[id-2]/dnu[id-1]);
        di1 = di[id-2]*::pow(anumin/dnu[id-2], a);
        
      }
      int ide = 22;
      if(anumax < dnu[BlzSimCommon::CSEED_SIZE-1]) {
        int idd;
        for(idd=id; idd<=BlzSimCommon::CSEED_SIZE; idd++) {
          if(anumax <= dnu[idd-1])
             break;
        }
        double a = log10(di[idd-2]/di[idd-1])/log10(dnu[idd-2]/dnu[idd-1]);
        die = di[idd-2]*::pow(anumax/dnu[idd-2], a);
        ide = idd;
      }
      else {
        die = di[BlzSimCommon::CSEED_SIZE-1];
      }
      double anui1 = anumin;
      double rat = anuf/(anui1*g1*g1);
      double ratr = 0.25 * rat;
      val1=(8.0+2.0*rat-rat*rat+4.0*rat*log(ratr))*(1.0e20/anui1)*(anuf/anui1)*di1*vala;
      if(val1 < 0.0)
        val1 = 0.0;
      //
      // Loop 1 to integrate over incoming photon frequency anui for lower gam value
      //
      for(int nu=id; nu<=ide; nu++) {
        double anui2 = dnu[nu-1];
        double di2 = di[nu-1];
        if(nu >= ide) {
          anui2=anumax;
          di2=die;
        }
        rat = anuf/(anui2*g1*g1);
        ratr = 0.25 * rat;
        val2=(8.0+2.0*rat-rat*rat+4.0*rat*log(ratr))*(1.0e20/anui2)*(anuf/anui2)*di2*vala;
        if(val2 < 0.0)
          val2=0.0;

        if((val1==0.0) || (val2==0.0)) {
          addit=0.5*(val1+val2)*(anui2-anui1);
          gran1 = gran1 + addit;
        }
        else {
          double test=abs((anui1-anui2)/anui1);
          if(test >= 0.001) {
            double ratnu = anui2/anui1;
            double a=1.0+log10(val2/val1)/log10(ratnu);
            double aa = abs(a);
            if((aa < 0.01) || (aa > 5.0)) {
              addit=0.5*(val1+val2)*(anui2-anui1);
            }
            else {
              addit=val1*(::pow(ratnu, a)-1.0) * anui1/a;
            }
            gran1 = gran1 + addit;
          }
        }

        anui1=anui2;
        val1=val2;
        di1=di2;
      } // End Loop 1 to integrate over incoming photon frequency anui for lower gam value
    } // if(anumin > anumax)

    double ratg = g2/g1;
    double ratgl = log10(ratg);
    valb = S0 * edist[ie]/(g2*g2);
    val1 = 0.0;
    val2 = 0.0;
    double gran2 = 0.0;
    //
    // Set up loop 2 to integrate over incident photon frequency anui
    //
    anumax=min(anuf,dnu[BlzSimCommon::CSEED_SIZE-1]);
    di1=di[0];
    id=2;
    anumin=0.25*anuf/(g1*g1);
    // skip over loop 2 entirely if this condition is not satisfied
    if(anumin <= anumax) {
      if(anumin <= dnu[0]) {
        anumin = dnu[0];
      }
      else {
        for(id=2; id<=BlzSimCommon::CSEED_SIZE; id++) {
          if(anumin <= dnu[id-1])
            break;
        }
        // if (anumin <= dnu[id-1]) is satisfied before id=CSEED_SIZE, then id will be set
        // to that number. Otherwise it will get to CSEED_SIZE (currently 22). The Fortran
        // code explicitly sets id to 22 in this case, but I don't see why - id will already
        // have the value 22
        double a = log10(di[id-2]/di[id-1])/log10(dnu[id-2]/dnu[id-1]);
        di1 = di[id-2]*::pow(anumin/dnu[id-2], a);
        
      }
      int ide = 22;
      if(anumax < dnu[BlzSimCommon::CSEED_SIZE-1]) {
        int idd;
        for(idd=id; idd<=BlzSimCommon::CSEED_SIZE; idd++) {
          if(anumax <= dnu[idd-1])
             break;
        }
        double a = log10(di[idd-2]/di[idd-1])/log10(dnu[idd-2]/dnu[idd-1]);
        die = di[idd-2]*::pow(anumax/dnu[idd-2], a);
        ide = idd;
      }
      else {
        die = di[BlzSimCommon::CSEED_SIZE-1];
      }
      double anui1 = anumin;
      double rat = anuf/(anui1*g2*g2);
      double ratr = 0.25 * rat;
      val1=(8.0+2.0*rat-rat*rat+4.0*rat*log(ratr))*(1.0e20/anui1)*(anuf/anui1)*di1*valb;
      if(val1 < 0.0)
        val1 = 0.0;
      //
      // Loop 2 to integrate over incoming photon frequency anui for lower gam value
      //
      for(int nu=id; nu<=ide; nu++) {
        double anui2 = dnu[nu-1];
        double di2 = di[nu-1];
        if(nu >= ide) {
          anui2=anumax;
          di2=die;
        }
        rat = anuf/(anui2*g2*g2);
        ratr = 0.25 * rat;
        val2=(8.0+2.0*rat-rat*rat+4.0*rat*log(ratr))*(1.0e20/anui2)*(anuf/anui2)*di2*valb;
        if(val2 < 0.0)
          val2=0.0;

        if((val1==0.0) || (val2==0.0)) {
          addit=0.5*(val1+val2)*(anui2-anui1);
          gran1 = gran1 + addit;
        }
        else {
          double test=abs((anui1-anui2)/anui1);
          if(test >= 0.001) {
            double ratnu = anui2/anui1;
            double a=1.0+log10(val2/val1)/log10(ratnu);
            double aa = abs(a);
            if((aa < 0.01) || (aa > 5.0)) {
              addit=0.5*(val1+val2)*(anui2-anui1);
            }
            else {
              addit=val1*(::pow(ratnu, a)-1.0) * anui1/a;
            }
            gran2 = gran2 + addit;
          }
        }

        anui1=anui2;
        val1=val2;
        di1=di2;
      } // End Loop 2 to integrate over incoming photon frequency anui for lower gam value
    } // if(anumin > anumax) Loop 2

    addit = 0.5 * (gran1 + gran2) * (g2 -g1);
    if((gran1!=0.0) && (gran2!=0.0)) {
      double a = 1.0 + log10(gran2/gran1)/ratgl;
      double aa = abs(a);
      if((aa>=.01) && (aa<=5.0)) {
        addit = gran1 * (::pow(ratg,a)-1.0)*g1/a;
      }
    }
    gran = gran + addit;
    g1 = g2;
    vala = valb;
  } // End loop over ie

  retVal = gran * 1e-20;
  return retVal;
 }

double BlzSim::ssc(const double anuf)
{
  // The Fortran is just a copy of ecdust() with some fairly minor modifications,
  // mainly the use of nuhi instead of BlzCommon::CDIST_SIZE. For expediency of porting, 
  // this C++ is likewise a copy of ecdust() above. At some point will consolidate
  // all this common code.
  double* gam = common.ggam;
  double* edist = common.edist;
  double* dnu = common.snu; // Note: this is different from ecdust()
  double* di = common.ssseed; // Note: this is different from ecdust()
  double retVal = 0.0, gran = 0.0;
  // (MSV: this is from the Fortran and the line of code is 
  // commented out in the Fortran, so including it all here for reference)
  // In this version, approximate that plasma velocity is along jet axis
  // and that Doppler factor of dust torus is the mean over its solid angle
  // as viewed in the plasma frame 
  // tdel=1.0/(gammad*(1.0d0-betad*dcsth1))
  double g1=gam[0];
  double vala=S0*edist[0]/(g1*g1);
  int ie;

  // prevent underflows
  if(di[common.nuhi-1] < 1.0e-25)
    common.nuhi--;

  //
  // Top-level loop to integrate over electron Lorentz factors (gam)
  //
  for(ie=1; ie<=BlzSimCommon::CDIST_SIZE-1; ie++) {
    double g2 = gam[ie], die;
    double valb = S0 * edist[ie]/(g2*g2);
    double val1=0.0, val2=0.0, gran1=0.0, addit=0.0;
    // Set up loop 1 to integrate over incident photon frequency anui
    double anumax=min(anuf,dnu[common.nuhi-1]);
    double di1=di[0];
    int id=2;
    double anumin=0.25*anuf/(g1*g1);
    // skip over loop 1 entirely if this condition is not satisfied
    if(anumin <= anumax) {
      if(anumin <= dnu[0]) {
        anumin = dnu[0];
      }
      else {
        for(id=2; id<=common.nuhi; id++) {
          if(anumin <= dnu[id-1])
            break;
        }
        // if (anumin <= dnu[id-1]) is satisfied before id=common.nuhi, then id will be set
        // to that number. Otherwise it will get to common.nuhi (e.g. 16). The Fortran
        // code explicitly sets id to common.nuhi in this case, but I don't see why - id will already
        // have the value common.nuhi
        di1 = 0.0;
        if((di[id-1]>=1.0e-25) && (di[id-2]>=1.0e-25)) {
          double a = log10(di[id-1]/di[id-2])/log10(dnu[id-1]/dnu[id-2]);
          di1 = di[id-1]*::pow(anumin/dnu[id-1], a);
        }
      }
      int ide = common.nuhi;
      if(anumax < dnu[common.nuhi-1]) {
        int idd;
        for(idd=id; idd<=common.nuhi; idd++) {
          if(anumax <= dnu[idd-1])
             break;
        }
        double a = log10(di[idd-2]/di[idd-1])/log10(dnu[idd-2]/dnu[idd-1]);
        die = di[idd-2]*::pow(anumax/dnu[idd-2], a);
        ide = idd;
      }
      else {
        die = di[common.nuhi-1];
      }
      double anui1 = anumin;
      double rat = anuf/(anui1*g1*g1);
      double ratr = 0.25 * rat;
      val1=(8.0+2.0*rat-rat*rat+4.0*rat*log(ratr))*(1.0e20/anui1)*(anuf/anui1)*di1*vala;
      if(val1 < 1.0e-28)
        val1 = 0.0;
      //
      // Loop 1 to integrate over incoming photon frequency anui for lower gam value
      //
      for(int nu=id; nu<=ide; nu++) {
        double anui2 = dnu[nu-1];
        double di2 = di[nu-1];
        if(nu >= ide) {
          anui2=anumax;
          di2=die;
        }
        rat = anuf/(anui2*g1*g1);
        ratr = 0.25 * rat;
        val2=(8.0+2.0*rat-rat*rat+4.0*rat*log(ratr))*(1.0e20/anui2)*(anuf/anui2)*di2*vala;
        if(val2 < 1.0e-28)
          val2=0.0;

        if((val1==0.0) || (val2==0.0)) {
          addit=0.5*(val1+val2)*(anui2-anui1);
          gran1 = gran1 + addit;
        }
        else {
          double test=abs((anui1-anui2)/anui1);
          if(test >= 0.001) {
            double ratnu = anui2/anui1;
            double a=1.0+log10(val2/val1)/log10(ratnu);
            double aa = abs(a);
            if((aa < 0.01) || (aa > 5.0)) {
              addit=0.5*(val1+val2)*(anui2-anui1);
            }
            else {
              addit=val1*(::pow(ratnu, a)-1.0) * anui1/a;
            }
            gran1 = gran1 + addit;
          }
        }

        anui1=anui2;
        val1=val2;
        di1=di2;
      } // End Loop 1 to integrate over incoming photon frequency anui for lower gam value
    } // if(anumin > anumax)

    double ratg = g2/g1;
    double ratgl = log10(ratg);
    valb = S0 * edist[ie]/(g2*g2);
    val1 = 0.0;
    val2 = 0.0;
    double gran2 = 0.0;
    //
    // Set up loop 2 to integrate over incident photon frequency anui
    //
    anumax=min(anuf,dnu[common.nuhi-1]);
    di1=di[0];
    id=2;
    anumin=0.25*anuf/(g1*g1);
    // skip over loop 2 entirely if this condition is not satisfied
    if(anumin <= anumax) {
      if(anumin <= dnu[0]) {
        anumin = dnu[0];
      }
      else {
        for(id=2; id<=common.nuhi; id++) {
          if(anumin <= dnu[id-1])
            break;
        }
        // if (anumin <= dnu[id-1]) is satisfied before id=common.nuhi, then id will be set
        // to that number. Otherwise it will get to common.nuhi (e.g. 16). The Fortran
        // code explicitly sets id to common.nuhi in this case, but I don't see why - id will already
        // have the value common.nuhi
        di1 = 0.0;
        addit = 0.0; // this doesn't occur in the 1st loop?
        if((di[id-1]>=1.0e-25) && (di[id-2]>=1.0e-25)) {
            double a = log10(di[id-1]/di[id-2])/log10(dnu[id-1]/dnu[id-2]);
            di1 = di[id-1]*::pow(anumin/dnu[id-1], a);
        }
      }
      int ide = common.nuhi;
      if(anumax < dnu[common.nuhi-1]) {
        int idd;
        for(idd=id; idd<=common.nuhi; idd++) {
          if(anumax <= dnu[idd-1])
             break;
        }
        double a = log10(di[idd-2]/di[idd-1])/log10(dnu[idd-2]/dnu[idd-1]);
        die = di[idd-2]*::pow(anumax/dnu[idd-2], a);
        ide = idd;
      }
      else {
        die = di[common.nuhi-1];
      }
      double anui1 = anumin;
      double rat = anuf/(anui1*g2*g2);
      double ratr = 0.25 * rat;
      val1=(8.0+2.0*rat-rat*rat+4.0*rat*log(ratr))*(1.0e20/anui1)*(anuf/anui1)*di1*valb;
      if(val1 < 1.0e-28)
        val1 = 0.0;
      //
      // Loop 2 to integrate over incoming photon frequency anui for lower gam value
      //
      for(int nu=id; nu<=ide; nu++) {
        double anui2 = dnu[nu-1];
        double di2 = di[nu-1];
        if(nu >= ide) {
          anui2=anumax;
          di2=die;
        }
        rat = anuf/(anui2*g2*g2);
        ratr = 0.25 * rat;
        val2=(8.0+2.0*rat-rat*rat+4.0*rat*log(ratr))*(1.0e20/anui2)*(anuf/anui2)*di2*valb;
        if(val2 < 1.0e-28)
          val2=0.0;

        if((val1==0.0) || (val2==0.0)) {
          addit=0.5*(val1+val2)*(anui2-anui1);
          gran1 = gran1 + addit;
        }
        else {
          double test=abs((anui1-anui2)/anui1);
          if(test >= 0.001) {
            double ratnu = anui2/anui1;
            double a=1.0+log10(val2/val1)/log10(ratnu);
            double aa = abs(a);
            if((aa < 0.01) || (aa > 5.0)) {
              addit=0.5*(val1+val2)*(anui2-anui1);
            }
            else {
              addit=val1*(::pow(ratnu, a)-1.0) * anui1/a;
            }
            gran2 = gran2 + addit;
          }
        }

        // The Fortran code does a goto to the next line *after* Loop 2, so a break should do the same here
        if(val2 < 1.0e-28)
          break;

        anui1=anui2;
        val1=val2;
        di1=di2;
      } // End Loop 2 to integrate over incoming photon frequency anui for lower gam value
    } // if(anumin > anumax) Loop 2

    addit = 0.5 * (gran1 + gran2) * (g2 -g1);
    if((gran1>=1.0e-28) && (gran2>=1.0e-28)) {
      double a = 1.0 + log10(gran2/gran1)/ratgl;
      double aa = abs(a);
      if((aa>=.01) && (aa<=5.0)) {
        addit = gran1 * (::pow(ratg,a)-1.0)*g1/a;
      }
    }
    gran = gran + addit;
    g1 = g2;
    vala = valb;
  } // End loop over ie

  double ssca = gran * 1e-20;
  if(ssca > 1.0e-30)
    retVal = ssca;

  return retVal;
 }

double BlzSim::polcalc(const double b, const double bx, const double by, const double bz, 
                       const double clos, const double slos)
{
  double chi = 0.0;
  double bxh=bx/b;
  double byh=by/b;
  double bzh=bz/b;
  double bdx = common.bdx;
  double bdy = common.bdy;
  double bdz = common.bdz;
  double gammad = common.gammad;
  double term1=slos*bxh+clos*bzh;
  double term2=slos*bdx+clos*bdz;
  double term3=(bdx*bxh+bdy*byh+bdz*bzh)*gammad/(1.0+gammad);
  double qx=bxh+(term1-term3)*bdx-term2*bxh;
  double qy=byh+(term1-term3)*bdy-term2*byh;
  double qz=bzh+(term1-term3)*bdz-term2*bzh;
  double q2=qx*qx+qy*qy+qz*qz;
  double ndq=qx*slos+qz*clos;
  term1 = -qy*clos;
  term2 = qx*clos-qz*slos;
  term3 = qy*slos;
  double term=sqrt(q2-ndq*ndq);
  double ex=term1/term; 
  double ey=term2/term;
  double ez=term3/term;
  term=ex*ex+ey*ey+ez*ez;
  double coschi=(-clos*ex+slos*ez)/sqrt(term);
  chi=acos(coschi);
  return chi;
}

void BlzSim::run(BlzSimInput& blzSimInput) 
{
  // TODO: fill this in :-)  
  // Seriously, this is where a most of the code ported from the "main" Fortran program will go, mostly as-is.
  // Then hopefully will have time to make it more modular after it is ported.
}
