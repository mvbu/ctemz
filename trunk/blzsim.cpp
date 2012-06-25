
#include <cmath>
#include <ctime>
#include "blzlog.h"
#include "blzmath.h"
#include "blzrand.h"
#include "blzsim.h"
#include "blzsiminputreader.h"

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
void BlzSim::psdsim(const int N, const double beta1, const double beta2, const double nu_break, const double t_incre1, double *lc_sim)
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
      flux_s[j-1]=::pow(nu[j-1], beta1);
      datareal[j-1]=sqrt((double)(fac_norm2*0.5*flux_s[j-1]))*R[j-1];
    }
    else {
      flux_s[j-1]=(::pow(nu[j-1], beta2)) * ((::pow(nu_break, beta1))/(::pow(nu_break, beta2)));
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
  double tdel=1.0/(pcm->gamd*(1.0D-pcm->betd*cs));
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
  // tdel=1.0/(gamd*(1.0d0-betad*dcsth1))
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
  // tdel=1.0/(gamd*(1.0d0-betad*dcsth1))
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
  double gammad = common.gamd;
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
  double ey=term2/term;
  if((ey > 1.0) && (ey < 1.000001))
    ey = .99999999;
  if((ey < -1.0) && (ey > -1.000001))
    ey = -.99999999;
  chi = asin(ey);
  return chi;
}

void BlzSim::initRandFromTime(bool bTestMode) 
{
  // Seeds the random number generator as closely as possible to what the Fortran code does
  time_t rawtime;
  struct tm* timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  // Fortran itime() gives 1-based values, C++ time() gives 0-based values, but both
  // seem to be returning the same numbers even though I'd expect the C++ numbers to
  // lower by 1. Since they don't seem to be, I will not add 1 to them.
  int hour = timeinfo->tm_hour; // + 1; // don't add 1 even though the docs would suggest that you should
  int min = timeinfo->tm_min; // + 1;
  int sec = timeinfo->tm_sec; // + 1;
  int ita = hour + min + sec; // the Fortran code just adds these up to create the seed
  randObj.setRandTestMode(bTestMode);
  double iseed = randObj.rand(ita);
  double dummy = randObj.rand(iseed);
}


void BlzSim::run(BlzSimInput& inp, double ndays, bool bTestMode)
{
  // This is where a most of the code ported from the "main" Fortran program will go, mostly as-is.
  // Then hopefully will have time to make it more modular after it is ported.
  const int D68=68;
  const int D100=100;
  const int D35=35;
  const int D1140=1140;
  const int D1141=1141;
  const int D44=44;
  const int D22=22;
  const int D4000=4000;
  const int D110=110;

  static int itarra[3];
  static double pq[D100][D1140][D68], pu[D100][D1140][D68], fpol[D100][D1140][D68],
    flux[D100][D1141][D68],gammin[D100][D1141],gammax[D100][D1141],
    xcell[D1141],gmax0[D100][D1141],egam[D100][D1141][44],
    phcell[D1141],rcell[D1141],ycell[D1141],cosph[D1141],
    sinph[D1141],bperp[D1141],bfield[D1141],n0[D100][D1141],
    zcell[D100][D1141],betadx[D100][D1141],betady[D100][D1141],
    betadz[D100][D1141],betad[D100][D1141],gammad[D100][D1141],
    betaux[D1141],betauy[D1141],betauz[D1141],
    betau[D1141],gammau[D1141],
    nu[D68],bx[D100][D1141],by[D100][D1141],bz[D100][D1141],
    delta[D100][D1141],enofe[D100][D1141][D44],fsync[D68],
    ididg[D1141],spsd[BLZSIM_DIM16384],
    gcnt[D44],igcnt[D44],fsynmd[D68][D4000],nouter[D1140],
    fsscmd[D68][D4000],fmdall[D68],deltmd[D1140],dmd[D1140],
    tlf[100],betamx[D1140],betamy[D1140],
    betamz[D1140],betamr[D1140],gamamr[D1140],bmdx[D4000],
    bmdy[D4000],bmdz[D4000],bmdtot[D4000],tlf1[D1140],
    flsync[D100][D1141][D68],flcomp[D100][D1141][D68],absorb[D68][D4000],
    fsync2[D68],cosmr[D1140],alphmd[D68][D4000],dustii[D110][D22],
    flec[D100][D1141][D68],flssc[D100][D1141][D68],mdd[D4000],useed[D110],
    phots[D68],phalph[D68],seedpk[D110],
    abexmd[D68][D4000],psi[D1140],sinpsi[D1140],cospsi[D1140],
    tanpsi[D1140];
  static float alfmdc[D68][D4000],syseed[D68],scseed[D68];
  int icelmx[D1141], imax[D1141], ibreak;
  double pol,pqcum,pucum,pmean,polc,ai2, pcum,tanv0,cosv0,gamup,beta,sinz,cosz,
    thlos,opang,tanop,cosop,sinop,zeta,tanz,slos,clos, eta,tanxi,xi,betacs,psiup,
    dth1,dth2,n0ave, betamd,betarl, cosmd,tanmd,cosbm,bup,bprp,ustob,tlfact,delt,
    gamb,glow,gmratl,gmratm,tloss,t1,t2,tlmin, eterm1,eterm2,glim,t2max,
    sinzps,coszps,tanzps,zetap, bd2,gm1,bmfact,n0mean;
  int dstart;

  // This line has no analog in the Fortran code, which has tdust 1) being filled in from the
  // input file (temzinp.txt) *and* 2) being part of a common block. Since I've created two
  // separate tdust's (one in BlzSimInput and one in BlzSimCommon), need to copy the value 
  // from BlzSimInput into BlzSimCommon. This need only be done once, and then the tdust in 
  // BlzSimCommon is what gets used/modified hereafter
  common.tdust = inp.tdust;

  // Initialize the random number generator. If bTestMode is true, the BlzRand::rand() method will get
  // numbers from a text file instead of actually generating new random numbers
  initRandFromTime(bTestMode);
  
  const int ICELLS = 50; // ICELLS cells along the axial direction
  const int NEND = 9; // NEND cells along each side of the hexagon
  int  jcells = 3*NEND*(NEND-1)+1; // jcells cells along the transverse direction (perp to axial dir)
  int  ancol = 2*NEND-1;
  double rbound = ancol*inp.rsize; // inp is the BlzSimInput object passed into this method
  int mdmax = 4000; // I think this matches the 4000 dimension in the arrays above
  // An SED will be printed out every ispecs time steps
  const int ISPECS = 1;
  int ispec = 1;
  // Set up frequencies
  int inu;
  for(inu=1; inu<=D68; inu++) {
    phots[inu-1] = 0.0;
    nu[inu-1] = ::pow(10.0, 10.0+0.25*(inu-1));
  }
  // alpha = (s-1)/2 is the underlying spectral index, where s = slope of electron E dist.
  // gmaxmn,gmaxmx are the min/max values of initial gamma_max of
  // electrons in a cell
  double gmaxmx = inp.gmrat * inp.gmaxmn;
  // gmin is the minimum value of initial gamma of electrons
  // 2p is the slope of the volume vs. initial gamma-max law over all cells,
  // V = V0*gamma_max**(-2p)
  double pexp = -2.0 * inp.p;
  double amppsd = 5.0;
  // Set up compilation of distribution of initial gamma_max values
  double gmrf=gmaxmx/inp.gmaxmn;
  double gmrfl=log10(gmrf)/43.0;
  int ie;
  for(ie=1; ie<=D44; ie++) {
    double gfl=log10(inp.gmaxmn)+gmrfl*(ie-1);
    gcnt[ie-1]=::pow(10.0, gfl);
    igcnt[ie-1] = 0;
  }
  const double FGEOM = 1.33;
  const double EMFACT = 3.74e-23; // not sure what this is, but code uses it (MSV June 2012)
  const double C6GAM = 8.2e-21;
  const double AMJY = 1.0e-26;
  common.zred1 = 1.0 + inp.zred;
  double sen = 2.0*inp.alpha + 1.0;
  zeta = inp.zeta/DEG_PER_RAD; // apparently inp.zeta is specified in DEG
  opang = inp.opang/DEG_PER_RAD; // ...and same with inp.opang
  // uratio is the ratio of energy density of electrons to that of mag. field
  // set the normalization of the electron energy distribution accordingly:
  n0ave = inp.uratio*inp.bave*inp.bave*(sen-2.0)/(8.0*PI*EMC2)/(::pow(inp.gmin, 2.0-sen) - ::pow(inp.gmaxmn,2.0-sen));
  // Line of sight
  thlos = inp.thlos/DEG_PER_RAD;
  slos = sin(thlos);
  clos = cos(thlos);
  sinz = sin(zeta);
  cosz = cos(zeta);
  tanz = sinz/cosz;
  gamup = 1.0/sqrt(1.0 - inp.betaup*inp.betaup);
  double betaus=inp.betaup*sinz;
  double gamus=1.0/sqrt(1.0-betaus*betaus); // I don't see where this is ever used in the Fortran code

  // Compression ratio of Mach disk. Ultra-relativistic eq. of state assumed, so compression ratio is that
  // given by Hughes, Aller, & Aller (1989, ApJ, 341, 54)
  double etac = sqrt(8.0*::pow(gamup, 4) - 17.0*gamup*gamup + 9.0)/gamup;
  double betaup2 = inp.betaup*inp.betaup;
  tanxi = (tanz*tanz*(3.0*betaup2-1.0) - (1-betaup2))/(tanz*(tanz*tanz + 1.0 + 2.0*betaup2));
  xi = atan(tanxi);
  betad[0][0] = sqrt(::pow(1.0 - ::pow(inp.betaup*cosz,2), 2) + 9.0*::pow(betaup2*cosz*sinz, 2))/(3.0*inp.betaup*sinz);
  gammad[0][0] = 1.0/sqrt(1.0-betad[0][0]*betad[0][0]);
  common.betd = betad[0][0];
  common.gamd = gammad[0][0];
  // Length of a cylindrical cell in pc
  double zsize = 2.0*inp.rsize/tanz;
  double volc = PI*inp.rsize*inp.rsize*zsize;
  // Length and volume of cell in plasma proper frame
  double zsizep= zsize/common.gamd;
  double volcp = volc/common.gamd;
  double svmd = sqrt(inp.vmd);
  double delobs = 1.0/(common.gamd*(1.0-common.betd*clos));
  dstart = mdmax+(ICELLS+100)*delobs+2500;
  // Time step in observer's frame in days
  double dtfact = (1.0-common.betd*clos)/(common.betd*clos);
  double dtime = 1190.0*zsize*dtfact*common.zred1;
  double itlast = ndays/dtime; // time "index" of the last timestep (quit when it is >= itlast)
  double mdrang = 0.5*(1.0/dtfact+1.0);
  // Distance of shock from axis and apex of conical jet
  tanop = tan(opang);
  cosop = cos(opang);
  sinop = tanop * cosop;
  // Next line is specific to the selected number of cells per slice
  double rshock = (2*NEND-1)*inp.rsize;
  // Distance of Mach disk from z value where conical shock intersects jet boundary
  double zshock = rshock/tanz;
  // Distance of Mach disk from from vertex of jet cone
  double zsvtex = rshock/tanop;
  double expon = inp.alpha+1.0;
  double exp1 = 0.5*expon-1.0;
  double anorm = (1.0+inp.alpha)/(inp.alpha+5.0/3.0);
  int id;
  // Computation of IR seed photon density from dust torus as function of distance down jet
  double zdist;

  for(id=1; id<=ICELLS+NEND; id++) {
    zdist = inp.zdist0+(id-NEND)*zsize+zshock;
    // Calculate min & max angles of dust torus in plasma frame
    double dphi1=asin(inp.dtrad/sqrt(zdist*zdist+inp.dtdist*inp.dtdist));
    double dphi2=atan(inp.dtdist/zdist);
    dth1 = dphi1 + dphi2;
    dth2 = dphi2 - dphi1;
    common.dcsth1 = -(cos(dth1) - common.betd)/(1.0-common.betd * cos(dth1));
    common.dcsth2 = -(cos(dth2) - common.betd)/(1.0-common.betd * cos(dth2));
    common.dsnth1 = sqrt(1.0 - common.dcsth1*common.dcsth1);
    common.dsnth2 = sqrt(1.0 - common.dcsth2*common.dcsth2);
    // Doppler factor of dust torus emission in frame of cells
    double tdel = 1.0/(common.gamd*(1.0-common.betd*common.dcsth1));
    common.dsang = TWOPI*(common.dsnth2 - common.dsnth1);
    // Calculate seed photon field from dust emission in plasma frame
    // Peak frequency of dust thermal emission for part of torus closest to shock
    seedpk[id-1] = 5.88e10*common.tdust*tdel;
    // Use this to set the frequency array of the seed photons
    for(inu=1; inu<=D22; inu++) {
      common.dustnu[inu-1] = seedpk[id-1]*::pow(10.0,-1.6+(inu-1)*0.1);
      common.dusti[inu-1] = seedph(common.dustnu[inu-1]);
      // If too far downstream, dust torus angles are in 2nd quadrant; need to correct
      if(common.dusti[inu-1] < 0.0)
        common.dusti[inu-1] = -common.dusti[inu-1];
      dustii[id-1][inu-1] = common.dusti[inu-1];
    }
  } // for(id=1; id<=ICELLS+NEND; id++)

  for(id=1; id<=ICELLS+NEND; id++) {
    double dflux=0.0;
    for(inu=2; inu<=D22; inu++) {
      if((dustii[id-1][inu-2]<=0.0) || (dustii[id-1][inu-1]<=0)) {
        dflux = dflux + 0.5*(dustii[id-1][inu-1]+dustii[id-1][inu-2]) * (common.dustnu[inu-1]-common.dustnu[inu-2]);
      }
      else {
        double a = log10(dustii[id-1][inu-1]/dustii[id-1][inu-2])/log10(common.dustnu[inu-1]/common.dustnu[inu-2]);
        dflux = dflux + dustii[id-1][inu-2]/(a+1.0)*common.dustnu[inu-2]*(::pow(common.dustnu[inu-1]/common.dustnu[inu-2], a+1.0) - 1.0);
      }

      zdist = inp.zdist0 + (id-NEND)*zsize+zshock;
    }
    useed[id-1] = dflux/C_CM_PER_SEC;
  }
  
  gmratl = log10(gmaxmx/inp.gmin)/40.0;

  double tinc = dtime;
  psdsim(BLZSIM_DIM16384, -inp.psdslp, -inp.psdslp, 1.0, tinc, spsd); 
  
  double psdsum = 0.0;
  int ip;
  for(ip=1; ip<=BLZSIM_DIM16384; ip++) {
    double spexp=1.0/(0.5*expon+1.0);
    spsd[ip-1] = ::pow(abs(spsd[ip-1]), spexp);
    psdsum = psdsum + (spsd[ip-1]/(double)BLZSIM_DIM16384);
  }

  for(ip=1; ip<=BLZSIM_DIM16384; ip++)
    spsd[ip-1]=amppsd*spsd[ip-1]/psdsum;
  
  int ip0 = randObj.rand(0) * 5000;
  int it = 0;
  // Set parameters of each cell at initial time
  // There are icells rows of cells, with jcells cells per row
  // The central cell is a Mach disk with low velocity; its observed radiation
  // is ignored, but it is an  important source of IC seed photons
  int i = 1, j = 0, nrow, ncol, neven, ncell, jold;
  for(nrow=-(NEND-1); nrow<=NEND-1; nrow++) {
    ncol = 2*NEND-(abs(nrow)+1);
    neven = ncol % 2;
    ncol = ncol/2;
    for(ncell=-ncol; ncell<=ncol; ncell++) {
      if((ncell==0) && (nrow==0)) {
        j=jcells;
        xcell[j-1]=0.0;
        ycell[j-1]=0.0;
        rcell[j-1]=0.0;
        cosph[j-1]=1.0;
        sinph[j-1]=0.0;
        j=jold;
      }
      else if((neven!=0) || (ncell!=0)) {
        j = j + 1;
        jold = j;
        xcell[j-1] = 2.0*ncell*inp.rsize;
        ycell[j-1] = nrow*SQRT3*inp.rsize;
        rcell[j-1] = sqrt(xcell[j-1]*xcell[j-1]+ycell[j-1]*ycell[j-1]);
        double zcol = rcell[j-1]/tanz;
        imax[j-1] = 2.0*zcol/zsize;
        if(imax[j-1] < 2)
          imax[j-1] = 2;
        // nouter[j-1] = approx. no. of cells between cell of interest and observer
        nouter[j-1]=imax[j-1];
        cosph[j-1]=xcell[j-1]/rcell[j-1];
        sinph[j-1]=ycell[j-1]/rcell[j-1];
        tanpsi[j-1]=rcell[j-1]/(zsvtex-zcol);
        psi[j-1] = atan(tanpsi[j-1]);
        cospsi[j-1] = cos(psi[j-1]);
        sinpsi[j-1]=tanpsi[j-1]*cospsi[j-1];
      }        
    } // for(ncell=-ncol; ncell<=ncol; ncell++)
  } //  for(nrow=(-NEND-1); nrow<=NEND-1)

  double zrf=zshock, xrf=0.0;

  // *** Set up Mach disk emission for earlier times ***
  // Compute time delay (no. of time steps) between when plasma passes Mach disk
  // and when it passes reference point of conical shock at column 30, in plasma frame
  double zmd = zrf - zshock;
  double idelmd = zmd/(zsize/common.betd);
  betadx[0][jcells-1] = 0.0;
  betady[0][jcells-1] = 0.0;
  betadz[0][jcells-1] = 1.0/3.0;
  betamd = betadz[0][jcells-1];
  double gammd = 1.0/sqrt(1.0 - betamd*betamd);
  betad[0][jcells-1] = betamd;
  gammad[0][jcells-1] = gammd;
  double dopref = common.gamd/gammd;
  double dopref2 = dopref*dopref;

  i = 1;
  int md;
  double xrand, phi;
  double gminmd = 0.15*inp.gmaxmn;

  // Determine B vector of MD assuming random magnetic field orientation
  for(md=1; md<=mdmax-1; md++) {
    xrand = randObj.rand(0);
    phi = TWOPI * xrand;
    xrand = randObj.rand(0);
    double costh = 2.0*(xrand-0.5);
    xrand = randObj.rand(0);
    double sign = xrand - 0.5;
    sign = sign/abs(sign);
    double thetab = sign*acos(costh);
    double sinthb=sin(thetab);
    double costhb=cos(thetab);
    double sinphb=sin(phi);
    double cosphb=cos(phi);
    // Compute B field components downstream of Mach disk shock
    int idelay = dstart-mdmax+md-idelmd;
    //if((idelay<1) || (idelay>(BLZSIM_DIM16384-ip0)))
    //  write(5,9222)idelay,i,j,md,ip0,zshock,zcell(i,j);
    double bavg = inp.bave*sqrt(spsd[idelay+ip0-1]);
    n0mean = n0ave*spsd[idelay+ip0-1];
    double n0prev = n0ave*spsd[idelay+ip0-2];
    j = jcells;
    n0[i-1][j-1] = etac*n0mean;
    n0prev = etac * n0prev;
    bx[i-1][j-1] = bavg*sinthb*cosphb*etac;
    by[i-1][j-1] = bavg*sinthb*sinphb*etac;
    bz[i-1][j-1] = bavg*costhb;
    bfield[j-1] = BlzMath::mag(bx[i-1][j-1], by[i-1][j-1], bz[i-1][j-1]);
    common.bfld = bfield[j-1];
    bmdx[md-1] = bx[i-1][j-1];
    bmdy[md-1] = by[i-1][j-1];
    bmdz[md-1] = bz[i-1][j-1];
    bmdtot[md-1] = common.bfld;
    // Calculate the initial maximum electron energy in the Mach disk
    bup = bavg*BlzMath::mag(costhb, gamup*sinthb);
    bprp = bavg*sqrt(::pow(costhb*sinz*sinz,2)+::pow(gamup*sinthb*cosz*cosz, 2) *
                     (::pow(cosphb*cosph[j-1]*cosph[j-1],2)+::pow(sinphb*sinph[j-1]*sinph[j-1], 2)));
    // Next line relates this to direction of B field relative to shock
    //  gmax0[i-1][j-1]=gmaxmn*(bprp/bup)**(2.0*pexp)
    // Next 3 lines assume a power-law distribution of gmax0, unrelated
    //     to direction of B field
    // See http://mathworld.wolfram.com/RandomNumber.html for choosing
    //  random numbers from a power-law distribution
    //  xrand=randproto(0)
    //  if(pexp.eq.0.0)gmax0[i-1][j-1]=gmaxmx
    //  if(pexp.lt.0.0)gmax0[i-1][j-1]=((gmaxmx**pexp-gmaxmn**pexp)*xrand+
    // ,  gmaxmn**pexp)**(1.0/pexp)
    gmax0[i-1][j-1] = 0.2*inp.gmaxmn;
    // Calculate energy distribution in the Mach disk cell
    // Compute energy density of photons in Mach disk, time delayed by 1 step
    // Ignore photons from dust torus since beaming is small in Mach disk
    int mdi = md-1;
    if(mdi <= 0)
      mdi = 1;
    // Value of SSC photon density for fast cooling case, from Sari and Esen (2001)
    double uphmd = bmdtot[mdi-1]*bmdtot[mdi-1]/(8.0*PI)*(sqrt(1.0+4.0*inp.uratio)-1.0)/2.0;
    ustob = 8.0*PI*uphmd/(bfield[j-1]*bfield[j-1]);
    tlfact = CC2*bfield[j-1]*bfield[j-1]*(1.0+ustob);
    // iend is the last slice of cells with energetic electrons
    int iend = i;
    delt = zsize*SEC_PER_YEAR*3.26/(gammd*betamd);
    gamb = gmax0[i-1][j-1]/(1.0+tlfact*gmax0[i-1][j-1]*delt);
    if(gamb < 1.0) gamb = 1.0;
    glow = gminmd/(1.0+tlfact*gminmd*delt);
    if(glow < 1.0) glow = 1.0;
    // Use 0.99 instead of 1 to avoid singularity of N(gamma) at gmax0
    inp.gmrat = 0.99*gmax0[i-1][j-1]/gminmd; // ick, Fortran is modifying the input variable
    gmratl = log10(inp.gmrat)/32.0; // TODO: copy inp.gmrat to a different variable, maybe common.gmrat
    gmratm = log10(gminmd/glow)/11.0;
    ibreak = 0;
    if(gamb <= glow) ibreak = 1;

    for(ie=1; ie<=D44; ie++) {
      egam[i-1][j-1][ie-1] = gminmd*::pow(10.0, gmratl*(ie-12));
      if(ie < 12)
        egam[i-1][j-1][ie-1] = glow*::pow(10.0, gmratm*(ie-1));
      if(ie == 1)
        egam[i-1][j-1][ie-1] = glow+0.2*(glow*::pow(10.0, gmratm)-glow);
      if((ibreak != 1) && (egam[i-1][j-1][ie-1] >= gamb)) {
        egam[i-1][j-1][ie-1] = gamb;
        ibreak = 1;
      }
      common.ggam[ie-1] = egam[i-1][j-1][ie-1];
      tloss = (gmax0[i-1][j-1] - common.ggam[ie-1])/(tlfact*common.ggam[ie-1]*gmax0[i-1][j-1]);
      t1 = (gminmd-common.ggam[ie-1])/(tlfact*common.ggam[ie-1]*gminmd);
      t2 = min(tloss,delt);
      tlmin = delt-(gminmd-common.ggam[ie-1])/(tlfact*gminmd*common.ggam[ie-1]);
      enofe[i-1][j-1][ie-1] = 0.0;
      eterm1 =1.0 - common.ggam[ie-1]*t1*tlfact;
      eterm2 =1.0 - common.ggam[ie-1]*t2*tlfact;
      if((eterm1>=0.0) && (eterm2>=0.0)) {
        enofe[i-1][j-1][ie-1] = n0[i-1][j-1]/((sen-1.0)*tlfact*::pow(common.ggam[ie-1],sen+1.0))
          * (1.0 - ::pow(eterm2, sen-1.0));
      }
      if(ie<12) {
        enofe[i-1][j-1][ie-1] = n0[i-1][j-1]/((sen-1.0)*tlfact*::pow(common.ggam[ie-1], sen+1.0))*
          (::pow(eterm1, sen-1.0) - ::pow(eterm2, sen-1.0));
      }
      //Divide by delt since integral is over time
      enofe[i-1][j-1][ie-1] = enofe[i-1][j-1][ie-1]/delt;
      if(enofe[i-1][j-1][ie-1]<0.0)
        enofe[i-1][j-1][ie-1] = 0.0;
      common.edist[ie-1] = enofe[i-1][j-1][ie-1];
    } // for(ie=1; ie<=D44; ie++)
    
    common.bperpp = common.bfld * (2.0/3.0);
    double restnu, ssabs, sstau, srcfn, fq1, fsyn1, absrb1;
    for(inu=1; inu<=40; inu++) { // TODO: what is the significance of 40??
      alfmdc[inu-1][md-1] = 10.0;
      fsscmd[inu-1][md-1] = 0.0;
      common.snu[inu-1] = nu[inu-1];
      restnu = nu[inu-1];
      // Synchrotron mean intensity for SSC calculation inside Mach disk
      ssabs = 1.02e4*(sen+2.0)*akapnu(restnu)*common.bperpp/(nu[inu-1]*nu[inu-1]);
      sstau = ssabs*CM_PER_PARSEC*inp.rsize*svmd;
      fsync[inu-1]=ajnu(restnu)*common.bperpp*inp.rsize*svmd*(EMFACT*CM_PER_PARSEC);
      if(fsync[inu-1] > 0.0) 
        common.nuhi=inu;
      if(sstau >= 0.01) {
        srcfn = fsync[inu-1]/sstau;
        if(sstau > 5.0)
          fsync[inu-1] = srcfn;
        if((sstau>0.01) && (sstau<=5.0))
          fsync[inu-1]= srcfn*(1.0-exp(-sstau));
      }
      // Now calculate synchrotron emission seen by other cells
      // Need to add a lower frequency to get spectral index of nu(1)
      if(inu <= 1) {
        fq1 = 0.98*nu[0];
        fsyn1 = ajnu(fq1/dopref)*dopref*dopref*common.bperpp;
        absrb1 = 1.02e4*(sen+2.0)*akapnu(fq1/dopref)*common.bperpp/::pow(fq1/dopref, 2);
      }
      fsynmd[inu-1][md-1] = ajnu(restnu/dopref)*dopref*dopref*common.bperpp;
      absorb[inu-1][md-1] = 1.02e4*(sen+2.0)*akapnu(restnu/dopref)*common.bperpp/(nu[inu-1]*nu[inu-1]);
      alphmd[inu-1][md-1] = 10.0;
      abexmd[inu-1][md-1] = 1.7;
      if((fsynmd[inu-1][md-1] > 0.0) && (fsyn1 > 0.0))
        alphmd[inu-1][md-1] = -log10(fsynmd[inu-1][md-1]/fsyn1)/log10(restnu/fq1);
      if((absorb[inu-1][md-1] > 0.0) && (absrb1 > 0.0))
        abexmd[inu-1][md-1] = -log10(fsynmd[inu-1][md-1]/fsyn1)/log10(restnu/fq1);
      // write(5,9994)md,inu,restnu,fsync[inu-1],(restnu/dopref),
      // ,   fsynmd[inu-1][md-1],sstau,
      // ,   absorb[inu-1][md-1],alphmd[inu-1][md-1],bperpp,dopref
      fq1 = restnu;
      fsyn1 = fsynmd[inu-1][md-1];
      common.ssseed[inu-1] = fsync[inu-1];
      // write(5,9911)md,j,dopref,bperpp,n0(i,j),ggam(44),edist(44),
      // ,  restnu,ssseed[inu-1],fsynmd[inu-1][md-1],fsscmd[inu-1][md-1]
    } // for(inu=1; inu<=40; inu++)

    for(inu=41; inu<=D68; inu++) {
      common.snu[inu-1] = nu[inu-1];
      alphmd[inu-1][md-1] = 10.0;
      common.ssseed[inu-1] = 0.0;
      fsync[inu-1] = 0.0;
      fsynmd[inu-1][md-1] = 0.0;
      absorb[inu-1][md-1] = 0.0;
    }

    // Calculate the SSC emission from the Mach disk for reference relative Doppler factor dopref
    fq1 = 0.98*nu[6];
    common.betd = betamd;
    common.gamd = gammd;
    double fssc1 = ssc(fq1/dopref)*dopref2;
    for(inu=7; inu<=D68; inu++) { // why 7?
      restnu = nu[inu-1];
      fsscmd[inu-1][md-1] = ssc(restnu/dopref)*dopref2;
      alfmdc[inu-1][md-1] = 10.0;
      if((fsscmd[inu-1][md-1] > 0.0) && (fssc1 > 0.0))
        alfmdc[inu-1][md-1] = -log10(fsscmd[inu-1][md-1]/fssc1)/log10(restnu/fq1);
      // write(5,9911)md,j,dopref,bfld,bperpp,n0(i,j),ggam(43),edist(43),
      // ,  restnu,ssseed(inu),fsynmd[inu-1][md-1],fsscmd[inu-1][md-1]
      fq1 = restnu;
      fssc1 = fsscmd[inu-1][md-1];
    }
  } // for(md=1; md<=mdmax-1; md++)

  //
  //     *** End Mach disk set-up ***
  // 

  i = 1;
  int ncells = 0;
  it = it + 1;
  md = mdmax;

  //
  // Start loop over all cells in first layer to set up physical parameters
  //
  for(j=1; j<=jcells; j++) {
    for(inu=1; inu<=D68; inu++) {
      fsynmd[inu-1][md-1] = 0.0;
      fsscmd[inu-1][md-1] = 0.0;
      fmdall[inu-1] = 0.0;
    }

    ididg[j-1] = 0;
    zcell[i-1][j-1] = zshock - (rcell[j-1]-inp.rsize)/tanz;
    int idelay = dstart + it - 1 + ((zrf-zcell[i-1][j-1]) +(xrf-xcell[j-1])*slos/(betad[0][0]*clos))/zsize;
    double bavg = inp.bave * sqrt(spsd[idelay+ip0-1]);
    n0mean = n0ave*spsd[idelay+ip0-1];
    // if(idelay.gt.(dstart+40).and.idelay.lt.(dstart+51));
    // ,  n0mean=10.0*n0mean
    //  write(5,9333)it,idelay,i,j,zlos,bavg,n0mean
    phcell[j-1] = atan2(sinph[j-1],cosph[j-1]);
    // Velocity vector of laminar component of pre-shock flow
    double betupx = inp.betaup*cosph[j-1]*sinpsi[j-1];
    double betupy = inp.betaup*sinph[j-1]*sinpsi[j-1];
    double betupz = inp.betaup*cospsi[j-1];
    // Velocity vector of the turbulent component of pre-shocked plasma
    phi = TWOPI*randObj.rand(0);
    double costh = 2.0*(randObj.rand(0)-0.5);
    double sign = randObj.rand(0)-0.5;
    sign = sign/abs(sign);
    double betaup2 = inp.betaup*inp.betaup;
    double thetat = sign*acos(costh);
    double sintht = sin(thetat);
    double costht = cos(thetat);
    double sinpht = sin(phi);
    double cospht = cos(phi);
    double betatx = inp.betat*cospht*sintht;
    double betaty = inp.betat*sinpht*sintht;
    double betatz = inp.betat*costht;
    double dotprd = betupx*betatx+betupy*betaty+betupz*betatz;
    double btparx = dotprd*betupx/betaup2;
    double btpary = dotprd*betupy/betaup2;
    double btparz = dotprd*betupz/betaup2;
    double btprpx = betatx-btparx;
    double btprpy = betaty-btpary;
    double btprpz = betatz-btparz;
    // Velocity vector of the pre-shock plasma including turbulent component
    double gamup2 = gamup*gamup;
    betaux[j-1] = (betupx+btparx+btprpx/gamup2)/(1.0+dotprd);
    betauy[j-1] = (betupy+btpary+btprpy/gamup2)/(1.0+dotprd);
    betauz[j-1] = (betupz+btparz+btprpz/gamup2)/(1.0+dotprd);
    betau[j-1] = BlzMath::mag(betaux[j-1], betauy[j-1], betauz[j-1]);
    gammau[j-1] = 1.0/sqrt(1.0-betau[j-1]*betau[j-1]);
    // Shock compression factor in terms of upstream speed and angle of shock zeta From Cawthorne & Cobb (1990)
    double sinzps = sin(zeta+psi[j-1]);
    double coszps = cos(zeta+psi[j-1]);
    double tanzps = sinzps/coszps;
    double eta = gammau[j-1]*betau[j-1]*sinzps*sqrt((8.0*::pow(betau[j-1]*sinzps, 2)-1.0/(gammau[j-1]*gammau[j-1]))/(1.0-::pow(betau[j-1]*coszps,2)));
    //  Deflection angle of velocity vector crossing shock, in rest frame of shock
    tanxi = (tanzps*tanzps*(3.0*betau[j-1]*betau[j-1]-1.0)-1.0/(gammau[j-1]*gammau[j-1]))/(tanzps*(tanzps*tanzps+1.0+2.0*betau[j-1]*betau[j-1]));
    xi = atan(tanxi);
    // Velocity vector downstream of shock
    double cszp2=coszps*coszps;
    betad[i-1][j-1] = sqrt(::pow(1.0-betau[j-1]*betau[j-1]*cszp2, 2)+9.0*::pow(betau[j-1],4)*cszp2*(1.0-cszp2))/(3.0*betau[j-1]*sqrt(1.0-cszp2));
    gammad[i-1][j-1] = 1.0/sqrt(1.0-betad[i-1][j-1]*betad[i-1][j-1]);
    double sinpsx = sin(psi[j-1]-xi);
    double cospsx = cos(psi[j-1]-xi);
    betadx[i-1][j-1] = betad[i-1][j-1]*cosph[j-1]*sinpsx;
    betady[i-1][j-1] = betad[i-1][j-1]*sinph[j-1]*sinpsx;
    betadz[i-1][j-1] = betad[i-1][j-1]*cospsx;
    // Determine B vector of cell assuming random magnetic field orientation
    // 11 continue (line from the Fortran code to help me orient myself :-)  )
    xrand = randObj.rand(0);
    phi = TWOPI * xrand;
    xrand = randObj.rand(0);
    costh = 2.0*(xrand-0.5);
    xrand = randObj.rand(0);
    sign = xrand-0.5;
    sign = sign/abs(sign);
    double thetab = sign*acos(costh);
    double sinthb = sin(thetab);
    double costhb = cos(thetab);
    double sinphb = sin(phi);
    double cosphb = cos(phi);

    // Compute B field components downstream of shock in the plasma frame
    // From equation derived by transforming from the upstream plasma
    // frame, compressing the component of B
    // parallel to the shock, and transforming the result to the rest frame
    // of the downstream plama
    double ssx, ssy, ssz, bux, buy, buz, bterm, bprpx, bprpy, bprpz;

    if(j != jcells) {
      zetap = atan(gammau[j-1]*tanzps)-atan(gammau[j-1]*tan(psi[j-1]));
      ssx = sin(zetap)*cosph[j-1];
      ssy = sin(zetap)*sinph[j-1];
      ssz = cos(zetap);
      bux = bavg*sinthb*cosphb;
      buy = bavg*sinthb*sinphb;
      buz = bavg*costhb;
      bterm = bux*ssx+buy*ssy+buz*ssz;
      // Finally, the downstream B field components
      bx[i-1][j-1] = bux+(eta-1.0)*bterm*ssx;
      by[i-1][j-1] = buy+(eta-1.0)*bterm*ssy;
      bz[i-1][j-1] = buz+(eta-1.0)*bterm*ssz;
      // Calculate upstream B field component perpendicular to the shock front bprp
      bup = sqrt(bux*bux+buy*buy+buz*buz);
      bprpx = bux-bterm*ssx;
      bprpy = buy-bterm*ssy;
      bprpz = buz-bterm*ssz;
      bprp = BlzMath::mag(bprpx, bprpy, bprpz);
    }
    else {
      // Set field of plasma in central cell, which is a Mach disk
      n0[i-1][j-1] = etac*n0mean;
      bx[i-1][j-1] = bavg*sinthb*cosphb*etac;
      by[i-1][j-1] = bavg*sinthb*sinphb*etac;
      bz[i-1][j-1] = bavg*costhb;
      bfield[j-1] = BlzMath::mag(bx[i-1][j-1], by[i-1][j-1], bz[i-1][j-1]);
      common.bfld = bfield[j-1];
      bmdx[md-1] = bx[i-1][j-1];
      bmdy[md-1] = by[i-1][j-1];
      bmdz[md-1] = bz[i-1][j-1];
      bmdtot[md-1] = common.bfld;
      // Calculate the initial maximum electron energy in the Mach disk
      bup = bavg*sqrt(costhb*costhb+::pow(gamup*sinthb, 2));
      bprp = bavg*sqrt(::pow(costhb*sinz*sinz,2)+::pow(gamup*sinthb*cosz*cosz,2)
                       *(::pow(cosphb*cosph[j-1]*cosph[j-1],2)+::pow(sinphb*sinph[j-1]*sinph[j-1],2)));
      //Next line relates this to direction of B field relative to shock
      //  gmax0[i-1][j-1]=gmaxmn*(bprp/bup)**(2.0*pexp)
      // Next 3 lines assume a power-law distribution of gmax0, unrelated
      //     to direction of B field
      // See http://mathworld.wolfram.com/RandomNumber.html for choosing
      //  random numbers from a power-law distribution
      xrand=randObj.rand(0);
      //  if(pexp.eq.0.0)gmax0[i-1][j-1]=gmaxmx
      //  if(pexp.lt.0.0)gmax0[i-1][j-1]=((gmaxmx**pexp-gmaxmn**pexp)*xrand+
      // ,  gmaxmn**pexp)**(1.0/pexp)
      gmax0[i-1][j-1] = 0.2*inp.gmaxmn;
      // gminmd = 0.15*inp.gmaxmn; gminmd already calculated earlier and value should not have changed by  ow

      // Calculate energy distribution in the Mach disk cell
      // Compute energy density of photons in Mach disk, time delayed by 1 step
      //   Ignore photons from dust torus since beaming is small in Mach disk
      int mdi = md-1;
      if(mdi <= 0)
        mdi = 1;
      double uphmd = bmdtot[mdi-1]*bmdtot[mdi-1]/(8.0*PI)*(sqrt(1.0+4.0*inp.uratio)-1.0)/2.0;
      ustob = 8.0*PI*uphmd/(bfield[j-1]*bfield[j-1]);
      tlfact = CC2*bfield[j-1]*bfield[j-1]*(1.0+ustob);
      // iend is the last slice of cells with energetic electrons
      int iend = i;
      delt = zsize*SEC_PER_YEAR*3.26/(gammd*betamd);
      gamb = gmax0[i-1][j-1]/(1.0+tlfact*gmax0[i-1][j-1]*delt);
      if(gamb < 1.0) gamb = 1.0;
      glow = gminmd/(1.0+tlfact*gminmd*delt);
      if(glow < 1.0) glow=1.0;
      inp.gmrat = 0.99*gmax0[i-1][j-1]/gminmd;
      gmratl = log10(inp.gmrat)/32.0;
      gmratm = log10(gminmd/glow)/11.0;
      ibreak = 0;
      if(gamb <= glow) ibreak = 1;
      //  write(5,9996)bfield[j-1],n0[i-1][j-1],gmax0[i-1][j-1],useed,usdmd,ustob,tlfact,
      // ,  delt,gamb,glow
      for(ie=1; ie<=D44; ie++) {
        egam[i-1][j-1][ie-1] = gminmd*::pow(10.0, gmratl*(ie-12));
        if(ie < 12)
          egam[i-1][j-1][ie-1] = glow*::pow(10.0, gmratm*(ie-1));
        if(ie == 1)
          egam[i-1][j-1][ie-1] = glow+0.2*(glow*::pow(10.0, gmratm)-glow);
        if((ibreak != 1) && (egam[i-1][j-1][ie-1] >= gamb)) {
          egam[i-1][j-1][ie-1] = gamb;
          ibreak = 1;
        }
        
        common.ggam[ie-1] = egam[i-1][j-1][ie-1];
        tloss = (gmax0[i-1][j-1]-common.ggam[ie-1])/(tlfact*common.ggam[ie-1]*gmax0[i-1][j-1]);
        t1=(gminmd-common.ggam[ie-1])/(tlfact*common.ggam[ie-1]*gminmd);
        t2 = min(tloss,delt);
        tlmin = delt-(gminmd-common.ggam[ie-1])/(tlfact*gminmd*common.ggam[ie-1]);
        enofe[i-1][j-1][ie-1] = 0.0;
        eterm1 = 1.0-common.ggam[ie-1]*t1*tlfact;
        eterm2 = 1.0-common.ggam[ie-1]*t2*tlfact;
        if((eterm1>=0.0) && (eterm2>=0.0)) {
          enofe[i-1][j-1][ie-1] = n0[i-1][j-1]/((sen-1.0)*tlfact*::pow(common.ggam[ie-1],sen+1.0))
            * (1.0 - ::pow(eterm2, sen-1.0));
        }
        if(ie<12) {
          enofe[i-1][j-1][ie-1] = n0[i-1][j-1]/((sen-1.0)*tlfact*::pow(common.ggam[ie-1], sen+1.0))*
            (::pow(eterm1, sen-1.0) - ::pow(eterm2, sen-1.0));
        }
        // Divide by delt since integral is over time
        enofe[i-1][j-1][ie-1] = enofe[i-1][j-1][ie-1]/delt;
        if(enofe[i-1][j-1][ie-1] < 0.0)
          enofe[i-1][j-1][ie-1] = 0.0;
        common.edist[ie-1] = enofe[i-1][j-1][ie-1];
      } // for(ie=1; ie<=D44; ie++)   1124 continue
    
      common.bperpp = common.bfld*(2.0/3.0);
      common.nuhi = 1;
      double restnu, ssabs, sstau, srcfn, fq1, fsyn1, absrb1;

      for(inu=1; inu<=40; inu++) {
        alfmdc[inu-1][md-1] = 10.0;
        fsscmd[inu-1][md-1] = 0.0;
        common.snu[inu-1] = nu[inu-1];
        restnu = nu[inu-1];
        // Synchrotron mean intensity for SSC calculation inside Mach disk
        ssabs = 1.02e4*(sen+2.0)*akapnu(restnu)*common.bperpp/(nu[inu-1]*nu[inu-1]);
        sstau = ssabs*CM_PER_PARSEC*inp.rsize*svmd;
        fsync[inu-1] = ajnu(restnu)*common.bperpp*inp.rsize*svmd*(EMFACT*CM_PER_PARSEC);
        if(fsync[inu-1] > 0.0)
          common.nuhi = inu;
        if(sstau >= 0.01) {
          srcfn = fsync[inu-1]/sstau;
          if(sstau > 5.0)
            fsync[inu-1] = srcfn;
          if((sstau>0.01) && (sstau<=5.0))
            fsync[inu-1] = srcfn*(1.0-exp(-sstau));
        } // 1126 continue

        // Now calculate synchrotron emission seen by other cells
        // Need to add a lower frequency to get spectral index of nu(1)
        if(inu <= 1) {
          fq1 = 0.98*nu[0];
          fsyn1 = ajnu(fq1/dopref)*dopref2*common.bperpp;
          absrb1 = 1.02e4*(sen+2.0)*akapnu(fq1/dopref)*common.bperpp/::pow(fq1/dopref,2);
        } //1127 continue
        
        fsynmd[inu-1][md-1] = ajnu(restnu/dopref)*dopref2*common.bperpp;
        absorb[inu-1][md-1] = 1.02e4*(sen+2.0)*akapnu(restnu/dopref)*common.bperpp/(nu[inu-1]*nu[inu-1]);
        alphmd[inu-1][md-1] = 10.0;
        abexmd[inu-1][md-1] = 1.7;
        if((fsynmd[inu-1][md-1]>0.0) && (fsyn1>0.0))
          alphmd[inu-1][md-1] = -log10(fsynmd[inu-1][md-1]/fsyn1)/log10(restnu/fq1);
        if((absorb[inu-1][md-1]>0.0) && (absrb1>0.0))
          abexmd[inu-1][md-1] = -log10(fsynmd[inu-1][md-1]/fsyn1)/log10(restnu/fq1);
        //  write(5,9996)restnu,fsync[inu-1],fsynmd[inu-1][md-1],sstau,
        // ,   absorb[inu-1][md-1],alphmd[inu-1][md-1],bperpp
        fq1 = restnu;
        fsyn1 = fsynmd[inu-1][md-1];
        common.ssseed[inu-1] = fsync[inu-1];
      } // 1125 continue

      for(inu=41; inu<=D68; inu++) { // do 1128 inu=41,68
        common.snu[inu-1] = nu[inu-1];
        alphmd[inu-1][md-1] = 10.0;
        common.ssseed[inu-1] = 0.0;
        fsync[inu-1] = 0.0;
        fsynmd[inu-1][md-1] = 0.0;
        absorb[inu-1][md-1] = 0.0;
      } // 1128 continue

      // Calculate the SSC emission from the Mach disk for reference
      //   relative Doppler factor dopref
      fq1 = 0.98*nu[6];
      common.betd = betamd;
      common.gamd = gammd;
      double fssc1 = ssc(fq1/dopref)*dopref2;

      for(inu=7; inu<=D68; inu++) { // do 1129 inu=7,68
        restnu = nu[inu-1];
        fsscmd[inu-1][md-1] = ssc(restnu/dopref)*dopref2;
        alfmdc[inu-1][md-1] = 10.0;
        if((fsscmd[inu-1][md-1]>0.0) && (fssc1>0.0))
          alfmdc[inu-1][md-1] = -log10(fsscmd[inu-1][md-1]/fssc1)/log10(restnu/fq1);
        //  write(5,9996)restnu,fsscmd[inu-1][md-1],alfmdc[inu-1][md-1]
        fq1 = restnu;
        fssc1 = fsscmd[inu-1][md-1];
      } // 1129 continue


      // Exit the for(j=1; j<=jcells; j++) loop
      break;
    }

    // 13 continue
    // Calculate component of magnetic field that is perpendicular to
    // the aberrated line of sight in the plasma frame
    // Line-of-sight vector in plasma frame
    bd2 = betad[i-1][j-1]*betad[i-1][j-1];
    gm1 = gammad[i-1][j-1]-1.0;
    double spx = (1.0+gm1*::pow(common.bdx/betad[i-1][j-1],2))*slos+gm1*common.bdx*common.bdz/bd2*clos;
    double spy = gm1*common.bdx*common.bdy/bd2*slos+gm1*common.bdy*common.bdz/bd2*clos;
    double spz = (1.0+gm1*::pow(common.bdz/betad[i-1][j-1],2))*clos+gm1*common.bdx*common.bdz/bd2*slos;
    //Need to normalize the line-of-sight vector - should be unit vector
    double spnorm = BlzMath::mag(spx, spy, spz);
    dotprd = (bx[i-1][j-1]*spx+by[i-1][j-1]*spy+bz[i-1][j-1]*spz)/spnorm;
    double bpx = bx[i-1][j-1]-dotprd*spx/spnorm;
    double bpy = by[i-1][j-1]-dotprd*spy/spnorm;
    double bpz = bz[i-1][j-1]-dotprd*spz/spnorm;
    bperp[j-1] = BlzMath::mag(bpx, bpy, bpz);
    bfield[j-1] = BlzMath::mag(bx[i-1][j-1], by[i-1][j-1], bz[i-1][j-1]);
    n0[i-1][j-1] = eta*n0mean;
    // Calculate the initial maximum electron energy in each cell
    // from the ratio of B(perp. to shock) to B_total in shock frame
    // Calculate upstream B field component perpendicular to the shock front bprp
    bprpx = bux-bterm*ssx;
    bprpy = buy-bterm*ssy;
    bprpz = buz-bterm*ssz;
    bprp = BlzMath::mag(bprpx, bprpy, bprpz);
    // Next line relates this to direction of B field relative to shock
    // gmax0[i-1][j-1]=gmaxmn*(bprp/bup)**(2.0*pexp)
    // Next 3 lines assume a power-law distribution of gmax0, unrelated
    // to direction of B field
    // see http://mathworld.wolfram.com/RandomNumber.html for choosing
    // random numbers from a power-law distribution
    xrand =randObj.rand(0);
    if(pexp == 0.0)
      gmax0[i-1][j-1] = gmaxmx;
    if(pexp<0.0)
      gmax0[i-1][j-1] = ::pow((::pow(gmaxmx,pexp)-::pow(inp.gmaxmn,pexp))*xrand+::pow(inp.gmaxmn,pexp), 1.0/pexp);
    
    int ig;
    for(ig=1; ig<=(D44-1); ig++) { // do 15 ig=1,43
      if((gmax0[i-1][j-1]<=gcnt[ig]) && (gmax0[i-1][j-1]>=gcnt[ig-1]))
        igcnt[ig-1] = igcnt[ig-1]+1;
    }      

  } // for(j=1; j<=jcells; j++)  End loop over all cells in FIRST LAYER   80 continue

  // Start another loop over all cells in first layer to calculate the
  // electron energy distribution and the emission
  double emisco, ecflux;

  for(j=1; j<=(jcells-1); j++) { // do 100 j=1,jcells-1
    ncells++;
    emisco = 0.0;
    ecflux = 0.0;
    zcell[i-1][j-1] = zshock-(rcell[j-1]-inp.rsize)/tanz;
    // Determine Doppler factor relative to the observer
    psiup = opang;
    common.bdx = betadx[i-1][j-1];
    common.bdy = betady[i-1][j-1];
    common.bdz = betadz[i-1][j-1];
    betacs = common.bdx*slos+common.bdz*clos;
    delta[i-1][j-1] = 1.0/(gammad[i-1][j-1]*(1.0-betacs));
    ecflux = 0.0;
    common.bperpp = bperp[j-1];
    common.bfld = bfield[j-1];
    // Determine velocity of the cell plasma relative to the MD plasma
    bmfact = 3.0*gammad[i-1][j-1]*betad[i-1][j-1]*betad[i-1][j-1];
    betamx[j-1] = (bmfact+(gammad[i-1][j-1]-1.0)*common.bdz)*common.bdx/(bmfact*(1.0+common.bdz/3.0));
    betamy[j-1] = (bmfact+(gammad[i-1][j-1]-1.0)*common.bdz)*common.bdy/(bmfact*(1.0+common.bdz/3.0));
    betamz[j-1] = ((bmfact+1.0)*betad[i-1][j-1]*betad[i-1][j-1]+(gammad[i-1][j-1]-1.0)*common.bdz*common.bdz)/(bmfact*(1.0+common.bdz/3.0));
    betamr[j-1] = BlzMath::mag(betamx[j-1], betamy[j-1], betamz[j-1]);
    gamamr[j-1] = 1.0/sqrt(1.0-(betamr[j-1]*betamr[j-1]));
    cosmr[j-1] = 0.0;
    double zcl = zcell[i-1][j-1];
    dmd[j-1] = sqrt(::pow(zshock-0.5*zsize-zcl, 2.0)+rcell[j-1]*rcell[j-1]);
    if(zshock != zcell[i-1][j-1]) {
      //  Determine Doppler factor of the MD emission in the cell's frame
      cosmr[j-1]=betamx[j-1]*xcell[j-1]+betamy[j-1]*ycell[j-1]+betamz[j-1]*(zshock-0.5*zsize-zcl)/gamamr[j-1];
      cosmr[j-1]=cosmr[j-1]/(betamr[j-1]*dmd[j-1]);
    }
    deltmd[j-1] = 1.0/(gamamr[j-1]*(1.0-betamr[j-1]*cosmr[j-1]));
    // Determine time step of Mach disk seed photons from light-travel delay
    double delcor = deltmd[j-1]/dopref;
    double mdmid = mdmax-(((zrf-zcl)*clos+(xrf-xcell[j-1])*slos)-dmd[j-1])/(dtfact*zsize);
    int md1 = mdmid-mdrang;
    if(md1 < 1) md1 = 1;
    int md2 = mdmid+mdrang;
    if(md2 > mdmax) md2 = mdmax;
    if(md1 > md2) md1 = md2;
    int amdrng = md2-md1+1.0;
    for(inu=1; inu<=D68; inu++) // do 3145 inu=1,68
      fmdall[inu-1] = 0.0;

    for(md=md1; md<=md2; md++) { // do 2147 md=md1,md2

      for(inu=1; inu<=D68; inu++) { 
        syseed[inu-1] = 0.0;
        scseed[inu-1] = 0.0;
        common.ssseed[inu-1] = 0.0;
      }

      double bp2cor = 0.0;

      // Cosine of Mach disk's B field in frame of cell plasma, from Lyutikov et al. (2003)
      double term1 = (bmdx[md-1]+bmdy[md-1]+bmdz[md-1])/(SQRT3*bmdtot[md-1]);
      double term2a = gamamr[j-1]/(bmdtot[md-1]*betamr[j-1])*(bmdx[md-1]*betamx[j-1]+bmdy[md-1]*betamy[j-1]+bmdz[md-1]*betamz[j-1]);
      double term2b = (gamamr[j-1]/(gamamr[j-1]+1.0))*(betamx[j-1]+betamy[j-1]+betamz[j-1])/SQRT3;
      double term3 = gamamr[j-1]*(1.0-(betamx[j-1]+betamy[j-1]+betamz[j-1])/SQRT3);
      cosbm = (term1+term2a*(term2b-1.0))/term3;
      // Correct for round-off error in case cosbm not within -1 to +1
      if(cosbm > 1.0) cosbm = 1.0;
      if(cosbm < -1.0) cosbm = -1.0;
      double bpcorr = sqrt(1.0-cosbm*cosbm)*1.5;

      // 
      // Calculate emission from various processes
      //

      // Calculate synchrotron flux from Mach disk in frame of cell plasma
      common.nuhi = 1;
      for(inu=1; inu<=40; inu++) { // do 146 inu=1,40
        if((fsynmd[inu-1][md-1]!=0.0) && (alphmd[inu-1][md-1]<=9.0)) {
          //  Synchrotron mean intensity for inverse Compton calculation
          double ssabs = absorb[inu-1][md-1]*::pow(bpcorr, 0.5*(sen+2.0))*::pow(delcor, abexmd[inu-1][md-1]);
          double sstau = ssabs*CM_PER_PARSEC*zsize;
          double ajofnu = fsynmd[inu-1][md-1]*::pow(bpcorr, 1.0+alphmd[inu-1][md-1])*::pow(delcor, 2.0+alphmd[inu-1][md-1]);
          ajofnu = ajofnu*zsize*(EMFACT*CM_PER_PARSEC);
          syseed[inu-1] = ajofnu;
          double srcfn = ajofnu/sstau;
          if(sstau > 5.0)
            syseed[inu-1] = srcfn;
          if((sstau>0.1) && (sstau<=5.0))
            syseed[inu-1] = srcfn*(1.0-exp(-sstau));
          syseed[inu-1] = syseed[inu-1]*volc*inp.vmd/(4.0*PI*zsize*dmd[j-1]*dmd[j-1]);
          double tauexp = sstau*dmd[j-1]/zsize;
          if(tauexp > 15.0)
            syseed[inu-1] = 0.0;
          if(tauexp <= 15.0)
            syseed[inu-1] = syseed[inu-1]/exp(tauexp);
        }

        scseed[inu-1] = 0.0;
        if(inu < 7)
          fmdall[inu-1] = fmdall[inu-1]+syseed[inu-1]/amdrng;
      } // for(inu=1; inu<=40; inu++)  146 continue

      //
      // Calculate inverse Compton flux from Mach disk in frame of cell plasma
      //
      for(inu=7; inu<=D68; inu++) { // do 147 inu=7,68
        // Inverse Compton mean intensity from Mach disk for 2nd-order
        // inverse Compton calculation in cells
        scseed[inu-1] = fsscmd[inu-1][md-1]*::pow(delcor, 2.0+alfmdc[inu-1][md-1]);
        scseed[inu-1] = scseed[inu-1]*volcp*inp.vmd/(4.0*PI*dmd[j-1]*dmd[j-1]);
        fmdall[inu-1] = fmdall[inu-1] + (syseed[inu-1]+scseed[inu-1])/amdrng;
      } // 147 continue

    } // for(md=md1; md<=md2; md++)  147 continue

    common.nuhi = 0;
    for(inu=1; inu<=D68; inu++) { // do 2148 inu=1,68
      common.ssseed[inu-1] = fmdall[inu-1];
      if(fmdall[inu-1] > 0.0)
        common.nuhi = inu;
    }

    // Calculate seed photon energy density energy loss calculation
    double usdmd = 0.0, aaa;
    for(inu=2; inu<=common.nuhi; inu++) { // do 149 inu=2,nuhi
      // Only calculate up to Klein-Nishina limit of gmin
      double epslon = 6.63e-27*common.snu[inu-1]/(inp.gmin*EMC2);
      if(epslon>1.0)
        break; // 150 continue
      if((fmdall[inu-1]>0.0) && (fmdall[inu-2]>0.0)) {
        aaa = log10(fmdall[inu-2]/fmdall[inu-1])/log10(common.snu[inu-1]/common.snu[inu-2]);
        usdmd = usdmd+4.0*PI/C_CM_PER_SEC/(1.0-aaa)*fmdall[inu-2]*common.snu[inu-2]*(::pow(common.snu[inu-1]/common.snu[inu-2],1.0-aaa)-1.0);
      }
    } //  149 continue

    // 150 continue
    // Calculate energy distribution in the cell
    id = (zcell[i-1][j-1]-zshock)/zsize+NEND;
    if(id < 1) id = 1;
    //if(id.gt.(icells+nend))write(3,9992)i,j,id,zcell[i-1][j-1],zshock,zsize
    ustob = 8.0*PI*(useed[id-1]+usdmd)/(bfield[j-1]*bfield[j-1]);
    tlfact =CC2*bfield[j-1]*bfield[j-1]*(1.0+ustob);
    tlf1[j-1] = tlfact;
    delt = zsize*SEC_PER_YEAR*3.26/(gammad[i-1][j-1]*betad[i-1][j-1]);
    gamb = gmax0[i-1][j-1]/(1.0+tlfact*gmax0[i-1][j-1]*delt);
    glow = inp.gmin/(1.0+tlfact*inp.gmin*delt);
    inp.gmrat = 0.99*gmax0[i-1][j-1]/inp.gmin;
    gmratl = log10(inp.gmrat)/32.0;
    gmratm = log10(inp.gmin/glow)/11.0;
    ibreak = 0;
    
    for(ie=1; ie<=D44; ie++) { //  do 90 ie=1,44
      egam[i-1][j-1][ie-1] = inp.gmin*::pow(10.0, gmratl*(ie-12));
      if(ie < 12)
        egam[i-1][j-1][ie-1] = glow*::pow(10.0,gmratm*(ie-1));
      if(ie == 1)
        egam[i-1][j-1][ie-1] = glow+0.2*(glow*::pow(10.0, gmratm)-glow);
      if((ibreak!=1) && (egam[i-1][j-1][ie-1]>=gamb)) {
        egam[i-1][j-1][ie-1] = gamb;
        ibreak = 1;
      }

      common.ggam[ie-1] = egam[i-1][j-1][ie-1];
      tloss = (gmax0[i-1][j-1]-common.ggam[ie-1])/(tlfact*common.ggam[ie-1]*gmax0[i-1][j-1]);
      t1 = (inp.gmin-common.ggam[ie-1])/(tlfact*common.ggam[ie-1]*inp.gmin);
      t2 = min(tloss,delt);
      tlmin = delt-(inp.gmin-common.ggam[ie-1])/(tlfact*inp.gmin*common.ggam[ie-1]);
      enofe[i-1][j-1][ie-1] = 0.0;
      eterm1 = 1.0-common.ggam[ie-1]*t1*tlfact;
      eterm2 = 1.0-common.ggam[ie-1]*t2*tlfact;
      if((eterm1>=0.0) && (eterm2>=0.0))
        enofe[i-1][j-1][ie-1] = n0[i-1][j-1]/((sen-1.0)*tlfact*::pow(common.ggam[ie-1], sen+1.0))*(1.0-::pow(eterm2,sen-1.0));
      if(ie< 12)
        enofe[i-1][j-1][ie-1] = n0[i-1][j-1]/((sen-1.0)*tlfact*::pow(common.ggam[ie-1], sen+1.0))*(::pow(eterm1,sen-1.0)-::pow(eterm2,sen-1.0));
      // Divide by delt since integral is over time
      enofe[i-1][j-1][ie-1] = enofe[i-1][j-1][ie-1]/delt;
      if(enofe[i-1][j-1][ie-1] < 0.0)
        enofe[i-1][j-1][ie-1] = 0.0;
      common.edist[ie-1] = enofe[i-1][j-1][ie-1];
    } // for(ie=1; ie<=D44; ie++)  90 continue

    delt = zsize*3.26*SEC_PER_YEAR/(gammad[i-1][j-1]*betad[i-1][j-1]);
    gammax[i-1][j-1] = gmax0[i-1][j-1]/(1.0+tlfact*delt*gmax0[i-1][j-1]);
    gammin[i-1][j-1] = inp.gmin/(1.0+tlfact*delt*inp.gmin);
    double emold = 0.0, chipol = 0.0;

    // Start loop over observed frequencies
    int ithin = 0;
    for(inu=1; inu<=D22; inu++) // do 93 inu=1,22
      common.dusti[inu-1] = dustii[id-1][inu-1];

    double spxec, spxssc, anumin, alnumn, emeold, emsold;
    int inumin, iinu;

    for(inu=1; inu<=D68; inu++)  { // do 95 inu=1,22
      double restnu = nu[inu-1]*common.zred1/delta[i-1][j-1];
      common.snu[inu-1] = nu[inu-1];
      double specin = 0.0001;
      emisco = 0.0;
      fsync2[inu-1] = 0.0;
      double ssabs = 0.0;
      ecflux = 0.0;
      double sscflx = 0.0;
      flsync[i-1][j-1][inu-1] = 0.0;
      flcomp[i-1][j-1][inu-1] = 0.0;
      flec[i-1][j-1][inu-1] = 0.0;
      flssc[i-1][j-1][inu-1] = 0.0;
      flux[i-1][j-1][inu-1] = 0.0;
      if(inu <= D44) {
        emisco = ajnu(restnu)*common.bperpp*delta[i-1][j-1]*delta[i-1][j-1];
        fsync2[inu-1] = emisco*zsize*(EMFACT*CM_PER_PARSEC);
        if(ithin != 1) {
          ssabs = 1.02e4*(sen+2.0)*akapnu(restnu)*common.bperpp/(nu[inu-1]*nu[inu-1])/delta[i-1][j-1];
          ssabs = ssabs*CM_PER_PARSEC*zsize*delta[i-1][j-1];
          //  Attenuate Mach disk emission from synchrotron absorption on way to cell
          double ssabsm = ssabs*dmd[j-1]/zsize;
          if(ssabsm >= 10.0)
            common.ssseed[inu-1] = 0.0;
          if((ssabsm<10.0) && (ssabsm > 0.1))
            common.ssseed[inu-1]= common.ssseed[inu-1]/exp(ssabsm);
          //  Return to absorption within cell
          // Use rsize instead of zsize because of aberration
          ssabs = ssabs*inp.rsize/zsize;
          if(ssabs<=(0.1/ancol))
            ithin = 1;
          double srcfn = fsync2[inu-1]/ssabs;
          if(ssabs>5.0)
            fsync2[inu-1] = srcfn;
          if((ssabs>0.1) && (ssabs<=5.0))
            fsync2[inu-1] = srcfn*(1.0-exp(-ssabs));
          double tauexp = nouter[j-1]*ssabs;
          if((rcell[j-1]>(0.98*rbound)) && (xcell[j-1]>=0.0))
            tauexp = 0.0;
          if(thlos==0.0)
            tauexp = nouter[j-1]*ssabs;
          // if(tauexp.gt.15.0)fsync2[inu-1] = 0.0;
          // if(tauexp.le.15.0)fsync2[inu-1]=fsync2[inu-1]/exp(tauexp)
          specin = inp.alpha;
        } // if(thin != 1)
        
        if(inu != 1) {
          if((emold>0.0) && (fsync2[inu-1]>0.0))
            specin = log10(emold/fsync2[inu-1])/log10(nu[inu-1]/nu[inu-2]);
        }

      } // if(inu <= D44) 92 continue

      flsync[i-1][j-1][inu-1] = fsync2[inu-1]*(volc/zsize)*common.zred1/(1.0e18*AMJY*inp.dgpc*inp.dgpc)*FGEOM;
      flux[i-1][j-1][inu-1] = flsync[i-1][j-1][inu-1];
      common.betd = betad[i-1][j-1];
      common.gamd = gammad[i-1][j-1];
      if(inu == 1) // looks like Fortran calculates polarization angle once - for the first frequency
        chipol = polcalc(bfield[j-1],bx[i-1][j-1],by[i-1][j-1],bz[i-1][j-1],clos,slos);
      if(specin < inp.alpha)
        specin = inp.alpha;
      double poldeg = (specin+1.0)/(specin+5.0/3.0);
      if(ssabs > 1.0)
        poldeg = 3.0/(12.0*specin+19);
      fpol[i-1][j-1][inu-1] = poldeg*flsync[i-1][j-1][inu-1];
      pq[i-1][j-1][inu-1] = fpol[i-1][j-1][inu-1]*cos(2.0*chipol);
      pu[i-1][j-1][inu-1] = fpol[i-1][j-1][inu-1]*sin(2.0*chipol);
      if(restnu >= 1.0e14) {
        spxec = 0.0001;
        spxssc = 0.0001;
        common.betd = betad[0][jcells-1];
        common.gamd = gammad[0][jcells-1];
        ecflux = ecdust(restnu)*3.086*volc*common.zred1*delta[i-1][j-1]*delta[i-1][j-1]/(inp.dgpc*inp.dgpc)*FGEOM;
        sscflx = ssc(restnu)*3.086*volc*common.zred1*delta[i-1][j-1]*delta[i-1][j-1]/(inp.dgpc*inp.dgpc)*FGEOM;
        double taupp = 0.0;
        if(nu[inu-1] >= 1.0e22) { // go to 99
          // Pair production opacity calculation
          // Expression for anumin includes typical interaction angle
          anumin = (1.24e20/(nu[inu-1]*common.zred1))*1.24e20*::pow(2.0*gammad[i-1][j-1], 2);
          alnumn = log10(anumin);
          inumin = ::pow(alnumn-10.0,4)+1;
          if(inumin <= 40) {  // go to 99
            if(anumin > nu[inumin-1])
              anumin = nu[inumin-1];

            for(iinu=inumin; iinu<=40; iinu++) { // do 97 iinu=inumin,40;
              double bp = sqrt(1.0-anumin/nu[iinu-1]);
              if(bp <= 0.0)
                bp=0.001;
              double xsecp1 = 1.25e-25*(1.0-bp*bp)*((3.0-::pow(bp,4))*log((1.0+bp)/(1.0-bp))-2.0*bp*(2.0-bp*bp));
              bp = sqrt(1.0-anumin/nu[iinu]);
              if(bp <= 0.0)
                bp=0.001;
              double xsecp2 = 1.25e-25*(1.0-bp*bp)*((3.0-::pow(bp,4))*log((1.0+bp)/(1.0-bp))-2.0*bp*(2.0-bp*bp));
              taupp = taupp+0.5*(phots[iinu-1]*xsecp1+phots[iinu]*xsecp2)*(nu[iinu]-nu[iinu-1])*nouter[j-1]*zsize*CM_PER_PARSEC;
            } // 97 continue

            if(taupp >= 0.1) { // go to 99
              if(taupp >= 10.0) {
                ecflux = ecflux/exp(taupp);
                sscflx = sscflx/exp(taupp);
              }
              else {
                ecflux = 0.0;
                sscflx = 0.0;
              }
            }
          } // if(inumin <= 40)
        } // if(nu[inu-1] >= 1.0e22)

        // 99
        flec[i-1][j-1][inu-1] = ecflux;
        flssc[i-1][j-1][inu-1] = sscflx;
        flcomp[i-1][j-1][inu-1] = ecflux+sscflx;
        flux[i-1][j-1][inu-1] = flsync[i-1][j-1][inu-1]+flcomp[i-1][j-1][inu-1];
        if((emeold>0.0) && (ecflux>0.0))
          spxec = log10(emeold/ecflux)/log10(nu[inu-1]/nu[inu-2]);
        if((emsold>0.0) && (sscflx>0.0))
          spxssc = log10(emsold/sscflx)/log10(nu[inu-1]/nu[inu-2]);
      } // if(restnu >= 1.0e14) 94 continue

      emold = fsync2[inu-1];
      emeold = ecflux;
      emsold = sscflx;
    } // for(inu=1; inu<=D68; inu++)  95 continue

    // 96 continue Fortran has this line here, but there's no reference to it (!)
    icelmx[j-1] = 1;

  } // for(j=1; j<=(jcells-1); j++) 100 continue

  // 
  // End loop over first layer of cells
  //

  if(icells != 1) {  // if(icells.eq.1)go to 2200
    i = i+1;
    istart = i;
  } // if(icells.eq.1)go to 2200   2200 continue
} // end run() method
