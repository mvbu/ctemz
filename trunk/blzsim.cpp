#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <string>
#include <omp.h>
#include "blzutil.h"
#include "blzlog.h"
#include "blzmath.h"
#include "blzrand.h"
#include "blzsim.h"
#include "blzsiminputreader.h"
#include "blztimer.h"

static const char FORMAT1[] = "%4d i %5d j %5d inu %5d ecflux %12.4E delta %12.4E\n";
static const char FORMAT2[] = "j %5d imax %5d\n";
static const char FORMAT3_1[] = "md%6d inu%5d fsync(inu)%15.3f\n";
static const char FORMAT3_2[] = "md%6d inu%5d fsynmd(inu,md)%15.2f\n";
static const char FORMAT3_3[] = "md%6d inu%5d fsscmd(inu,md)%12.5E\n";
static const char FORMAT14001[] = "%5d i%5d j%5d md%6d inu%5d %s%10.4f\n";
static const char FORMAT14005[] = "%5d i%5d j%5d md%6d inu%5d %s%10.4f %s%10.4f\n";
static const char FORMAT_ARRAY_INT[] = "%7d %8d\n"; // format 14006 in temz.f
static const char FORMAT_ARRAY_FLOAT[] = "%7d %12.5E\n"; // format 14007 in temz.f
static const char FORMAT_EDIST[] = "%5s %7d %7d %2d %8d %12.5E %12.5E %12.3E %12.3E %12.4E %12.5E\n"; // format 14008 in temz.f
static const char FORMAT_STRING_FLOAT[] = "%10s %8d %12.5E\n"; // format 14009 in temz.f
static const char FORMAT_STRING_3INT_FLOAT6_FLOAT8[] = "%10s %7d %7d %7d %12.6E %12.8E\n";

const double SMALL_FMDALL = 1e-30;

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

/* Output/result parameter lc_sim must be of dimension BLZSIM_DIM131072 */
void BlzSim::psdsim(const int N, const double beta1, const double beta2, const double nu_break, const double t_incre1, double *lc_sim)
{
  // This code is ported from Fortran routine temz.f:psdsim() by R. Chatterjee
  // The fortran code contains a bunch of stuff like dataim_s1[] and flux_s2[] that 
  // don't get used at all, so I did not include them here
  double nu[BLZSIM_DIM131072];
  //double dat[BLZSIM_DIM262144],R[BLZSIM_DIM262144];
  double *dat,*R;
  double dataim[BLZSIM_DIM131072],datareal[BLZSIM_DIM131072],flux_s[BLZSIM_DIM131072];
  //double nu[BLZSIM_DIM65536],dat[BLZSIM_DIM131072],R[BLZSIM_DIM131072];
  //double dataim[BLZSIM_DIM65536],datareal[BLZSIM_DIM65536],flux_s[BLZSIM_DIM65536];
  double fac_norm,fac_norm2;
  float  ann;
  int j,i,nn,ISEED1,ISEED2,isign;
  dat = (double*)malloc(BLZSIM_DIM262144*sizeof(double));
  R = (double*)malloc(BLZSIM_DIM262144*sizeof(double));

  fac_norm=1./(N*t_incre1);
  fac_norm2=pow((double)N,2)/(2.*N*t_incre1);

  // randObj is an instance of BlzRand inside this BlzSim instance
  ISEED1=393521;
  ISEED2=17263;
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

  free(R);

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

  free(dat);
}

double BlzSim::scatcs(const double vx, const double vy, const double vz,
                      const double sx, const double sy, const double sz, 
                      const double x, const double y, const double z)
{
  double g,g1,g2,slx,sly,slz,smx,smy,smz;
  double s, cscat=0.0;
  // v is velocity in units of c, s is line-of-sight unit vector
  // (x,y,z) is position of cell relative to Mach disk, cscat is
  // cosine of scattering angle in plasma frame
  double v2 = BlzMath::magSquared(vx, vy, vz);
  g2 = 1.0/(1.0-v2);
  g = sqrt(g2);
  g1 = g-1.0;
  // Determine unit vector of l.o.s. in plasma frame (slx,sly,slz)
  slx=(1.0+g1*vx*vx/v2)*sx+g1*vx*vy*sy/v2+g1*vx*vz*sz/v2;
  sly=g1*vx*vy*sx/v2+(1.0+g1*vy*vy/v2)*sy+g1*vy*vz*sz/v2;
  slz=g1*vx*vz*sx/v2+g1*vy*vz*sy/v2+(1.0+g1*vz*vz/v2)*sz;
  s = BlzMath::mag(slx, sly, slz);
  slx=slx/s;
  sly=sly/s;
  slz=slz/s;
  // Determine unit vector from Mach disk to cell in plasma frame (smx,smy,smz)
  smx=(1.0+g1*vx*vx/v2)*x+g1*vx*vy*y/v2+g1*vx*vz*z/v2;
  smy=g1*vx*vy*x/v2+(1.0+g1*vy*vy/v2)*y+g1*vy*vz*z/v2;
  smz=g1*vx*vz*x/v2+g1*vy*vz*y/v2+(1.0+g1*vz*vz/v2)*z;
  s = BlzMath::mag(smx, smy, smz);
  smx=smx/s;
  smy=smy/s;
  smz=smz/s;
  // Calculate cosine of scattering angle in plasma frame using head-on
  // approximation (photon is scattering in direction of electron motion)
  cscat = slx*smx+sly*smy+slz*smz;
  return cscat;
}

double BlzSim::xseckn(const double epsi, const double epsf, const double g, const double y)
{
  // Subroutine to calculate the Klein-Nishina cross-section; from Dermer & Menon
  // (2009, High Energy Radiation from Black Holes: Gamma Rays, Cosmic Rays,
  // and Neutrinos, Princeton U. Press), eq. 6.31
  double z = epsf/(g*epsi*y);
  double xsec = (y+1.0/y-2.0*z+z*z)/(g*epsi);
  return xsec;
}

double BlzSim::seedph(const double f)
{
  double cs1, cs2;
  cs1 = common.csth2;
  cs2 = common.csth1;
  common.freq = f;
  double csm = 0.5*(cs1+cs2);
  double ans1 = BlzMath::qg5(cs1, csm, sdgran, &common);
  double ans2 = BlzMath::qg5(csm, cs2, sdgran, &common);
  // Multiply by 2pi, but then divide by 4pi to get intensity
  double retVal = ans1 + ans2;
  return retVal;
}

// Non-member function
double sdgran(double cs, void *pObject)
{
  BlzSimCommon* pcm = (BlzSimCommon*)pObject;
  const double HOK = 4.80e-11;
  double retVal = 0.0;
  double csp = -(cs-pcm->betd)/(1.0-pcm->betd*cs);
  double tdel=1.0/(pcm->gamd*(1.0-pcm->betd*csp));
  double f = pcm->freq/tdel;
  double expon= HOK * f/pcm->tdust;
  double val = 0.0;

  // The Fortran code has a bunch of goto statements. The structure
  // below implements the same logic, I think :-)
  if(expon <= .01) {
    val=(2.76e-16*f)*pcm->tdust*(f/9.0e20);
  }
  else if(expon > 10.0) {
    if(expon > 20.0)
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
  double g1=0;
  int ie, ie1=1;
  double ef = HOMC2*anuf;
  double ef1 = ef + 1;
  double xsecc=0.0; // , xsecx=0.0; // xsecx used only in write statements that have not been ported here

  // Set the highest frequency of seed photons for scattering
  // This block is here in the same spot as in the Fortran, but seems like it could eventually be moved
  // to after all the ef and gam checks. TODO: check if this is possible and move it. -msv May2013
  int nhi = BlzSimCommon::CSEED_SIZE, nuhi = BlzSimCommon::CSEED_SIZE, nn;
  for(nn=10; nn<=nhi; nn++) {
    if((dnu[nn-1] >= anuf) || (di[nn-1]<=1.0e-30)) { //  go to 2010
      nuhi = nn - 1;
      break; // we have our value of nuhi: jump out of the loop
    }
  }

  if(ef1 <= gam[0]) {
    g1 = gam[0];
  }
  else if(ef1 >= gam[BlzSimCommon::CDIST_SIZE-1]) {
    return 0.0;
  }
  else if(gam[BlzSimCommon::CDIST_SIZE-1]/gam[0] < 1.01) {
    return 0.0;
  }
  else {
    for(ie=2; ie<=BlzSimCommon::CDIST_SIZE; ie++) {
      if(gam[ie-1] > ef1)
        break;
    }
    // Fortran sets ie=CDIST_SIZE here if the above "break" is never hit, but that should be unnecessary since
    // ie should already be CDIST_SIZE if the break is never hit. -msv Jan2013
    ie1 = ie-1;
    g1 = ef1;
  }

  do {
    if(gam[ie1]/gam[ie1-1] > 1.0002)
      break;
    ie1++;
    if(ie1 >= BlzSimCommon::CDIST_SIZE-1)
      return 0.0;
  } while(ie1 < BlzSimCommon::CDIST_SIZE-1);

  double S0 = S0_ECDUST;
  double vala = S0*edist[ie1-1]/g1; //13

  //
  // Top-level loop to integrate over electron Lorentz factors (gam)
  //
  for(ie=ie1; ie<=BlzSimCommon::CDIST_SIZE-1; ie++) { // do 3000 ie=ie1,43
    double g2 = gam[ie], die;
    double valb = S0 * edist[ie]/g2;
    double val1=0.0, val2=0.0, gran1=0.0, addit=0.0;
    // Set up loop 1 to integrate over incident photon frequency anui
    double anumax=min(anuf,dnu[nuhi-1]);
    double di1=di[0];
    int id=2;
    double anumin = anuf/((2.0*g1)*(g1-ef)*(1.0-common.csang*sqrt(1.0-1.0/(g1*g1))));

    // skip over loop 1 entirely if this condition is not satisfied
    if(anumin < anumax) { // go to 601
      if(anumin <= dnu[0]) {
        anumin = dnu[0];
      }
      else {
        for(id=2; id<=nuhi; id++) { // do 3
          if(anumin <= dnu[id-1])
            break;
        } // 3 continue
        // if (anumin <= dnu[id-1]) is satisfied before id=nuhi, then id will be set
        // to that number. Otherwise it will get to nuhi. The Fortran
        // code explicitly sets id to nuhi in this case, but I don't see why - id will already
        // have the value nuhi
        double a = log10(di[id-2]/di[id-1])/log10(dnu[id-2]/dnu[id-1]);
        di1 = di[id-2]*::pow(anumin/dnu[id-2], a);
        
      }
      int ide = nuhi;
      if(anumax < dnu[nuhi-1]) {
        int idd;
        for(idd=id; idd<=nuhi; idd++) {
          if(anumax <= dnu[idd-1])
             break;
        }
        double a = log10(di[idd-2]/di[idd-1])/log10(dnu[idd-2]/dnu[idd-1]);
        die = di[idd-2]*::pow(anumax/dnu[idd-2], a);
        ide = idd;
      }
      else {
        die = di[nuhi-1];
      }
      double anui1 = anumin;
      double xx = 1.0-ef/g1;
      double rat = HOMC2*anui1*g1*(1.0-common.csang*sqrt(1.0-1.0/(g1*g1)));
      xsecc = xseckn(rat,ef,g1,xx);
      xsecc = xsecc*(g1*HOMC2*anui1);
      val1 = xsecc*(1.0e20/anui1)*(anuf/anui1)*di1*vala;

      if(val1 < 1.0e-40)
        val1 = 0.0;

      // 25 continue (don't think the Fortran ever references this "25" label)

      //
      // Loop 1 to integrate over incoming photon frequency anui for lower gam value
      //
      for(int nu=id; nu<=ide; nu++) { // do 600 nu=id,ide
        double anui2 = dnu[nu-1];
        double di2 = di[nu-1];
        if(nu >= ide) {  // go to 30
          anui2=anumax;
          di2=die;
          xx = 1.0-ef/g1;
        }

        // 30
        rat = HOMC2*anui2*g1*(1.0-common.csang*sqrt(1.0-1.0/(g1*g1)));
        xsecc = xseckn(rat,ef,g1,xx);
        xsecc = xsecc*(g1*HOMC2*anui2);
        val2 = xsecc*(1.0e20/anui2)*(anuf/anui2)*di2*vala;

        if(val2 < 1.0e-40)
          val2 = 0.0;

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
      } // End Loop 1 to integrate over incoming photon frequency anui for lower gam value 600 continue
    } // if(anumin > anumax)  601 continue

    double ratg = g2/g1;
    double ratgl = log10(ratg);
    valb = S0 * edist[ie]/g2;
    val1 = 0.0;
    val2 = 0.0;
    double xx;
    double gran2 = 0.0;

    //
    // Set up loop 2 to integrate over incident photon frequency anui
    //
    anumax=min(anuf,dnu[nuhi-1]);
    di1=di[0];
    id=2;
    anumin = anuf/((2.0*g2)*(g2-ef)*(1.0-common.csang*sqrt(1.0-1.0/(g2*g2))));
    // skip over loop 2 entirely if this condition is not satisfied
    if(anumin < anumax) {
      if(anumin <= dnu[0]) {
        anumin = dnu[0];
      }
      else {
        for(id=2; id<=nuhi; id++) { // do 1003
          if(anumin <= dnu[id-1])
            break;
        }
        // if (anumin <= dnu[id-1]) is satisfied before id=nuhi, then id will be set
        // to that number. Otherwise it will get to nuhi. The Fortran
        // code explicitly sets id to nuhi in this case, but I don't see why - id will already
        // have the value nuhi
        di1 = 0.0;
        addit = 0.0;
        if((di[id-1]>=1.0e-25) && (di[id-2]>=1.0e-25)) { // go to 1005
          double a = log10(di[id-1]/di[id-2])/log10(dnu[id-1]/dnu[id-2]);
          di1 = di[id-1]*::pow(anumin/dnu[id-1], a);
        }
        
      }

      int ide = nuhi; // 1005
      if(anumax < dnu[nuhi-1]) {
        int idd;
        for(idd=id; idd<=nuhi; idd++) {
          if(anumax <= dnu[idd-1])
             break;
        }
        double a = log10(di[idd-2]/di[idd-1])/log10(dnu[idd-2]/dnu[idd-1]);
        die = di[idd-2]*::pow(anumax/dnu[idd-2], a);
        ide = idd;
      }
      else {
        die = di[nuhi-1];
      }

      // 1009 continue
      double anui1 = anumin;
      xx = 1.0-ef/g2;
      double rat = HOMC2*anui1*g2*(1.0-common.csang*sqrt(1.0-1.0/(g2*g2)));
      xsecc = xseckn(rat,ef,g2,xx);
      xsecc = xsecc*(g2*HOMC2*anui1);
      val1 = xsecc*(1.0e20/anui1)*(anuf/anui1)*di1*valb;

      if(val1 < 1.0e-40)
        val1 = 0.0;

      gran2 = 0.0; // 1025

      //
      // Loop 2 to integrate over incoming photon frequency anui for lower gam value
      //
      if(id < ide) { // skip over loop 2 if this condition is not met
        for(int nu=id; nu<=ide; nu++) { // do 1600 nu=id,ide
          double anui2 = dnu[nu-1];
          double di2 = di[nu-1];
          if(nu >= ide) {
            anui2=anumax;
            di2=die;
          }

          // 1030
          rat = HOMC2*anui2*g2*(1.0-common.csang*sqrt(1.0-1.0/(g2*g2)));
          xsecc = xseckn(rat,ef,g2,xx);
          xsecc = xsecc*(g2*HOMC2*anui2);
          val2 = xsecc*(1.0e20/anui2)*(anuf/anui2)*di2*valb;

          if(val2 < 1.0e-40)
            val2=0.0;

          if((val1<1.0e-38) || (val2<1.0e-38)) {
            addit=0.5*(val1+val2)*(anui2-anui1);
            gran2 = gran2 + addit;
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
        }  // 1600 continue. End anui loop 2 (End Loop 2 to integrate over incoming photon frequency anui for lower gam value)
      } // if(id < ide)
    } // if(anumin > anumax) Loop 2

    // 1601 continue
    addit = 0.5 * (gran1 + gran2) * (g2 -g1);
    if((gran1>1e-28) && (gran2>1e-28) && (ratgl>=0.0001)) {
      double a = 1.0 + log10(gran2/gran1)/ratgl;
      double aa = abs(a);
      if((aa>=.01) && (aa<=5.0)) {
        addit = gran1 * (::pow(ratg,a)-1.0)*g1/a;
      }
    }
    gran = gran + addit;
    g1 = g2;
    vala = valb;
  } // 3000 continue. End gam loop (End loop over ie)

  retVal = gran * 1e-20;
  return retVal;
 }

double BlzSim::ssc(const double anuf)
{
  // The Fortran is just a copy of ecdust() with some fairly minor modifications,
  // mainly the use of nuhi instead of BlzCommon::CDIST_SIZE. For expediency of porting, 
  // this C++ is likewise a copy of ecdust() above. At some point will consolidate
  // all this common code.
  double* dnu = common.snu; // Note: this is different from ecdust()
  double* di = common.ssseed; // Note: this is different from ecdust()
  double retVal = 0.0, gran = 0.0;
  double xsecc=0.0; // , xsecx=0.0; // xsecx used only in write statements that have not been ported here

  // Set the highest frequency of seed photons for scattering
  // This block is here in the same spot as in the Fortran, but seems like it could eventually be moved
  // to after all the ef and gam checks. TODO: check if this is possible and move it. -msv May2013
  int nhi = BlzSimCommon::CSSC_SIZE, nn;
  int nuhi = BlzSimCommon::CSSC_SIZE;

  for(nn=20; nn<=nhi; nn++) {
    if((common.snu[nn-1] >= anuf) || (common.ssseed[nn-1]<=1.0e-30)) {
      nuhi = nn - 1;
      break; // we have our value of nuhi: jump out of the loop
    }
  }

  double g1=common.ggam[0];
  double vala;
  int ie, ie1=1;
  double ef = HOMC2*anuf;
  double ef1 = ef + 1.0;

  if(ef1 <= common.ggam[0]) {
    g1 = common.ggam[0];
  }
  else if(ef1 >= common.ggam[BlzSimCommon::CDIST_SIZE-1]) {
    return 0.0;
  }
  else if(common.ggam[BlzSimCommon::CDIST_SIZE-1]/common.ggam[0] < 1.01) {
    return 0.0;
  }
  else {
    for(ie=2; ie<=BlzSimCommon::CDIST_SIZE; ie++) {
      if(common.ggam[ie-1] > ef1)
        break;
    }
    // Fortran sets ie=CDIST_SIZE here if the above "break" is never hit, but that should be unnecessary since
    // ie should already be CDIST_SIZE if the break is never hit. -msv Jan2013
    ie1 = ie-1;
    g1 = ef1;
  }

  do {
    if(common.ggam[ie1]/common.ggam[ie1-1] > 1.0002)
      break;
    ie1++;
    if(ie1 >= BlzSimCommon::CDIST_SIZE-1)
      return 0.0;
  } while(ie1 < BlzSimCommon::CDIST_SIZE-1);
  
  double S0 = S0_SSC;
  vala = S0*common.edist[ie1-1]/g1;

  //
  // Top-level loop to integrate over electron Lorentz factors (gam)
  //
  for(ie=ie1; ie<=BlzSimCommon::CDIST_SIZE-1; ie++) {
    double g2 = common.ggam[ie], die;
    double valb = S0 * common.edist[ie]/g2;
    double val1=0.0, val2=0.0, gran1=0.0, addit=0.0;
    // Set up loop 1 to integrate over incident photon frequency anui
    double anumax=min(anuf,dnu[nuhi-1]);
    double di1=di[0];
    int id=2;
    double anumin = anuf/((2.0*g1)*(g1-ef)*(1.0-common.cscat*sqrt(1.0-1.0/(g1*g1))));
    // skip over loop 1 entirely if this condition is not satisfied
    if(anumin < anumax) { // go to 601
      if(anumin <= dnu[0]) {
        anumin = dnu[0];
      }
      else {
        for(id=2; id<=nuhi; id++) {
          if(anumin <= dnu[id-1])
            break;
        }
        // if (anumin <= dnu[id-1]) is satisfied before id=nuhi, then id will be set
        // to that number. Otherwise it will get to nuhi (e.g. 16). The Fortran
        // code explicitly sets id to nuhi in this case, but I don't see why - id will already
        // have the value nuhi
        di1 = 0.0;
        if((di[id-1]>=1.0e-25) && (di[id-2]>=1.0e-25)) {
          double a = log10(di[id-1]/di[id-2])/log10(dnu[id-1]/dnu[id-2]);
          di1 = di[id-1]*::pow(anumin/dnu[id-1], a);
        }
      }
      int ide = nuhi;
      if(anumax < dnu[nuhi-1]) {
        int idd;
        for(idd=id; idd<=nuhi; idd++) {
          if(anumax <= dnu[idd-1])
             break;
        }
        double a = log10(di[idd-2]/di[idd-1])/log10(dnu[idd-2]/dnu[idd-1]);
        die = di[idd-2]*::pow(anumax/dnu[idd-2], a);
        ide = idd;
      }
      else {
        die = di[nuhi-1];
      }

      // 9 continue
      double anui1 = anumin;
      double xx = 1.0-ef/g1;

      double rat = HOMC2*anui1*g1*(1.0-common.cscat*sqrt(1.0-1.0/(g1*g1)));
      xsecc = xseckn(rat,ef,g1,xx);
      xsecc = xsecc*(g1*HOMC2*anui1);
      val1 = xsecc*(1.0e20/anui1)*(anuf/anui1)*di1*vala;

      if(val1 < 1.0e-40)
        val1 = 0.0;

      // 25 continue (don't think the Fortran ever references this "25" label)

      //
      // Loop 1 to integrate over incoming photon frequency anui for lower gam value
      //
      for(int nu=id; nu<=ide; nu++) { // do 600 nu=id,ide
        double anui2 = dnu[nu-1];
        double di2 = di[nu-1];
        if(nu >= ide) { // go to 30
          anui2=anumax;
        }

        // 30

        rat = HOMC2*anui2*g1*(1.0-common.cscat*sqrt(1.0-1.0/(g1*g1)));
        xsecc = xseckn(rat,ef,g1,xx);
        xsecc = xsecc*(g1*HOMC2*anui2);
        val2 = xsecc*(1.0e20/anui2)*(anuf/anui2)*di2*vala;

        if(val2 < 1.0e-40)
          val2 = 0.0;

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

    // 601 continue
    double ratg = g2/g1;
    double ratgl = log10(ratg);
    valb = S0 * common.edist[ie]/g2;
    val1 = 0.0;
    val2 = 0.0;
    double gran2 = 0.0;
    double xx;

    //
    // Set up loop 2 to integrate over incident photon frequency anui
    //
    anumax=min(anuf,dnu[nuhi-1]);
    di1=di[0];
    id=2;
    anumin = anuf/((2.0*g2)*(g2-ef)*(1.0-common.cscat*sqrt(1.0-1.0/(g2*g2))));

    // skip over loop 2 entirely if this condition is not satisfied
    if(anumin <= anumax) {
      if(anumin <= dnu[0]) {
        anumin = dnu[0];
      }
      else {
        for(id=2; id<=nuhi; id++) {
          if(anumin <= dnu[id-1])
            break;
        }
        // if (anumin <= dnu[id-1]) is satisfied before id=nuhi, then id will be set
        // to that number. Otherwise it will get to nuhi (e.g. 16). The Fortran
        // code explicitly sets id to nuhi in this case, but I don't see why - id will already
        // have the value nuhi
        di1 = 0.0;
        addit = 0.0; // this doesn't occur in the 1st loop?
        if((di[id-1]>=1.0e-25) && (di[id-2]>=1.0e-25)) {
            double a = log10(di[id-1]/di[id-2])/log10(dnu[id-1]/dnu[id-2]);
            di1 = di[id-1]*::pow(anumin/dnu[id-1], a);
        }
      }
      int ide = nuhi;
      if(anumax < dnu[nuhi-1]) {
        int idd;
        for(idd=id; idd<=nuhi; idd++) {
          if(anumax <= dnu[idd-1])
             break;
        }
        double a = log10(di[idd-2]/di[idd-1])/log10(dnu[idd-2]/dnu[idd-1]);
        die = di[idd-2]*::pow(anumax/dnu[idd-2], a);
        ide = idd;
      }
      else {
        die = di[nuhi-1];
      }
      double anui1 = anumin;
      xx = 1.0-ef/g2;
      double rat = HOMC2*anui1*g2*(1.0-common.cscat*sqrt(1.0-1.0/(g2*g2)));
      xsecc = xseckn(rat,ef,g2,xx);
      xsecc = xsecc*(g2*HOMC2*anui1);
      val1 = xsecc*(1.0e20/anui1)*(anuf/anui1)*di1*valb;

      if(val1 < 1.0e-40)
        val1 = 0.0;

      gran2 = 0.0;

      //
      // Loop 2 to integrate over incoming photon frequency anui for lower gam value
      //
      if(id < ide) {
        for(int nu=id; nu<=ide; nu++) { // do 1600 nu=id,ide
          double anui2 = dnu[nu-1];
          double di2 = di[nu-1];
          if(nu >= ide) {
            anui2=anumax;
            di2=die;
          }

          // 1030
          rat = HOMC2*anui2*g2*(1.0-common.cscat*sqrt(1.0-1.0/(g2*g2)));
          xsecc = xseckn(rat,ef,g2,xx);
          xsecc = xsecc*(g2*HOMC2*anui2);
          val2 = xsecc*(1.0e20/anui2)*(anuf/anui2)*di2*valb;

          if(val2 < 1.0e-40)
            val2=0.0;

          if((val1==0.0) || (val2==0.0)) {
            addit=0.5*(val1+val2)*(anui2-anui1);
            gran2 = gran2 + addit;
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

          // can't find the equiv of these lines in the latest temz.f (msv 1/11/2013), so commenting out
          // The Fortran code does a goto to the next line *after* Loop 2, so a break should do the same here
          //if(val2 < 1.0e-28)
          //  break;

          anui1=anui2;
          val1=val2;
          di1=di2;
        } // 1600 continue. End anui loop2 (End Loop 2 to integrate over incoming photon frequency anui for lower gam value)
      } // if(id < ide)
    } // if(anumin > anumax) Loop 2

    // 1601
    addit = 0.5 * (gran1 + gran2) * (g2 -g1);
    if((gran1>=1.0e-28) && (gran2>=1.0e-28) && (ratgl>=0.0001)) {
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

//
//    vdcalc computes downstream velocity vector; Lorentz transformation of unit
//    vector along shock front follows Lyutikov et al. (2003, ApJ, 597, 998)
//
void BlzSim::vdcalc(const double vx, const double vy,const double vz,
                    const double sx, const double sy,const double sz,
                    double* vdx, double* vdy, double* vdz, double *vd, // outputs
                    double* gd, double* eta)
{
  double  v2,v,g,spx,spy,spz,s,g2;
  double  dotprd,vparx,vpary,vparz,vprpx,vprpy,vprpz,vprp2,vprp,uprp;
  double  vd2,vfact, vdprpx,vdprpy,vdprpz,vdprp2,gdprp;
  vfact = 1.0;
  v2 = BlzMath::magSquared(vx, vy ,vz);
  v = sqrt(v2);
  g2 =1.0/(1.0-v2);
  g = sqrt(g2);
  dotprd = vx*sx+vy*sy+vz*sz;
  // Calculate components of velocity vector parallel and perpendicular to shock front
  vparx=dotprd*sx;
  vpary=dotprd*sy;
  vparz=dotprd*sz;
  vprpx=vx-vparx;
  vprpy=vy-vpary;
  vprpz=vz-vparz;
  vprp2=vprpx*vprpx+vprpy*vprpy+vprpz*vprpz;
  vprp=sqrt(vprp2);
  // Shock jump condition, from Konigl (1980, Phys. Fluids, 23, 1083)
  uprp=g*vprp;
  // uprp must exceed the proper sound speed, 1/sqrt(2) for a shock
  // Otherwise, it is a sound wave and the velocity does not change significantly
  if(uprp > 0.7071)
    vfact=(1.0+1.0/(g2*vprp2))/(3.0);
  vdprpx=vprpx*vfact;
  vdprpy=vprpy*vfact;
  vdprpz=vprpz*vfact;
  vdprp2=BlzMath::magSquared(vdprpx, vdprpy, vdprpz);
  // Shock compression ratio (downstream density/upstream density)
  // From Cawthorne + Cobb (1990, ApJ, 350, 536)
  *eta = sqrt((8.0*v2*common.sinz*common.sinz-1.0/g2)/(1.0-v2*common.cosz*common.cosz))*v*g*common.sinz;
  // eta must be >= 1; if formula gives < 1, it is not a shock
  if(*eta < 1.0)
    *eta = 1.0;
  *vdx=vparx+vdprpx;
  *vdy=vpary+vdprpy;
  *vdz=vparz+vdprpz;
  vd2=BlzMath::magSquared(*vdx, *vdy, *vdz);
  *vd=sqrt(vd2);
  *gd=1.0/sqrt(1.0-vd2);
}

//
//    bdcalc() computes magnetic field components downstream of shock
//
void BlzSim::bdcalc(const double vx, const double vy, const double vz, 
                    const double ax, const double ay, const double az, // unit vector normal to shock
                    const double bx, const double by, const double bz, // magnetic field vector
                    double* v, double* g,  // outputs
                    double* bparx, double* bpary, double* bparz,
                    double* bprpx, double* bprpy, double* bprpz,
                    double* bpar, double* bprp)
{
  double dotprd,sx,sy,sz,a,v2,g1,g2,gp1,apx,apy,apz;
  // v is velocity in units of c, a is a unit vector normal to the shock front
  // b is magnetic field vector
  v2 = BlzMath::magSquared(vx, vy, vz);
  *v = sqrt(v2);
  g2 = 1.0/(1.0-v2);
  *g = sqrt(g2);
  g1 = *g-1.0;
  gp1 = *g+1.0;
  // Determine unit vector of shock normal. in plasma frame
  dotprd = vx*ax+vy*ay+vz*az;
  apx = *g*ax-g1*dotprd*vx/v2;
  apy = *g*ay-g1*dotprd*vy/v2;
  apz = *g*az-g1*dotprd*vz/v2;
  a = BlzMath::mag(apx, apy, apz);
  apx = apx/a;
  apy = apy/a;
  apz = apz/a;
  // Compute component of B that is normal to the shock front
  dotprd = bx*apx+by*apy+bz*apz;
  // Calculate components of B field parallel and perpendicular to shock normal
  *bprpx = dotprd*apx;
  *bprpy = dotprd*apy;
  *bprpz = dotprd*apz;
  *bprp = BlzMath::mag(*bprpx, *bprpy, *bprpz);
  *bparx = bx-*bprpx;
  *bpary = by-*bprpy;
  *bparz = bz-*bprpz;
  *bpar = BlzMath::mag(*bparx, *bpary, *bparz);
}

//
//    bcalc computes magnetic field component parallel and perpendicular to shock
//    front or line of sight; follows Lyutikov et al. (2003, ApJ, 597, 998)
//
void BlzSim::bcalc(const double vx, const double vy, const double vz,
                   const double sx, const double sy, const double sz,
                   const double bx, const double by, const double bz,
                   const double v, const double g,
                   double* bparx, double* bpary, double* bparz, // outputs
                   double* bprpx, double* bprpy, double* bprpz,
                   double* bpar, double* bprp)
{
  double s,v2,g1,g2,gp1,denom,spx,spy,spz;
  double b,dotprd;
  // v is velocity in units of c, s is line-of-sight or shock front unit vector
  // b is magnetic field vector
  v2 = v*v;
  g2 = g*g;
  g1 = g-1.0;
  gp1 = g+1.0;
  // Determine unit vector of shock front or l.o.s. in plasma frame
  dotprd = vx*sx+vy*sy+vz*sz;
  spx = g*sx-g1*dotprd*vx/v2;
  spy = g*sy-g1*dotprd*vy/v2;
  spz = g*sz-g1*dotprd*vz/v2;
  s = BlzMath::mag(spx, spy, spz);
  spx = spx/s;
  spy = spy/s;
  spz = spz/s;
  // Calculate components of B field parallel and perpendicular to shock front
  // or l.o.s.
  dotprd = bx*spx+by*spy+bz*spz;
  *bparx = dotprd*spx;
  *bpary = dotprd*spy;
  *bparz = dotprd*spz;
  *bpar = BlzMath::mag(*bparx, *bpary, *bparz);
  *bprpx = bx-*bparx;
  *bprpy = by-*bpary;
  *bprpz = bz-*bparz;
  *bprp = BlzMath::mag(*bprpx, *bprpy, *bprpz);
}

double adjustTrig(double in) 
{
  double retVal = in;
  if(in >= 1.0) retVal=0.999999;
  else if(in <= -1.0) retVal=-0.999999;
  return retVal;
}

int BlzSim::getIwp(int it)
{
  int iwp = 0;
  if(it==1) iwp = 1;
  if(it==25 || it==50) iwp = 1;
  if(it==75 || it==100) iwp = 1;
  if(it==125 || it==150) iwp = 1;
  if(it==175 || it==200) iwp = 1;
  if(it==225 || it==250) iwp = 1;
  if(it==275 || it==300) iwp = 1;
  if(it==325 || it==350) iwp = 1;
  if(it==375 || it==400) iwp = 1;
  if(it==425 || it==450) iwp = 1;
  if(it==475 || it==500) iwp = 1;
  if(it==525 || it==550) iwp = 1;
  if(it==575 || it==600) iwp = 1;
  if(it==625 || it==650) iwp = 1;
  if(it==675 || it==700) iwp = 1;
  if(it==725 || it==750) iwp = 1;
  if(it==775 || it==800) iwp = 1;
  if(it==825 || it==850) iwp = 1;
  if(it==875 || it==900) iwp = 1;
  if(it==925 || it==950) iwp = 1;
  if(it==975 || it==1000) iwp = 1;

  return iwp;
}

void writeedist(FILE* fp, double* arr, int nsize)
{
  fprintf(fp, "%10s\n", "edist");
  int i;
  for(i=0; i<nsize; i++) {
    fprintf(fp, "%12.5f\n",arr[i]);
  }
}

// Max nTestOut is 8. Default is 0 (i.e. no test output)
void BlzSim::run(BlzSimInput& inp, double ndays, bool bTestMode, bool bSingleThreaded, int nTestOut)
{
  int nThreads = bSingleThreaded ? 1 : BlzUtil::getNumProcessors();
  const int NUM_THREADS = nThreads;
  BlzLog::warnScalar("BlzSim::run() NUM_THREADS = ", NUM_THREADS);
  BlzLog::warnScalar("BlzSim::run() nTestOut = ", nTestOut);
  // Create object to see how much time ssc() takes to execute
  BlzTimer sscTimer(NUM_THREADS);
  BlzTimer parallelTimer(1);
  // This is where a most of the code ported from the "main" Fortran program will go, mostly as-is.
  // Then hopefully will have time to make it more modular after it is ported.
  const int D68=68;
  const int D400=400;
  const int D35=35;
  const int D1140=1140;
  const int D1141=1141;
  const int D44=44;
  const int D22=22;
  const int D130000=130000;
  const int D130000PADDED=131000;
  const int D451=451;

  // calculate and store starting points for each thread in the loops that deal with frequencies 1 to 68, or 7 to 68
  int threadIntervals7[NUM_THREADS][2];
  BlzMath::getSubIntervals(7, D68, NUM_THREADS, threadIntervals7);
  int threadIntervals1[NUM_THREADS][2];
  BlzMath::getSubIntervals(1, D68, NUM_THREADS, threadIntervals1);
  BlzIndexTracker indexTracker7(7, D68);
  BlzIndexTracker indexTracker1(1, D68);

  bool bOutputFilesCreated = false;
  FILE* pfSpec = NULL; // 3 ctemzspec.txt
  FILE* pfLc = NULL; // 4 ctemzlc.txt
  FILE* pfPol = NULL; // 5 ctemzpol.txt
  FILE* pfMap = NULL; // 7 ctemzpol.txt
  FILE* pfTestOut = NULL; // 9 ctestout.txt
  std::string specFile("ctemzspec.txt");
  std::string lcFile("ctemzlc.txt");
  std::string mapFile("ctemzmap.txt"); // this will change once "it" value is known/set
  std::string polFile("ctemzpol.txt");
  std::string testOutFile("ctestoutp.txt");
  std::string dpath("maps/");
  if(nTestOut>0) pfTestOut = fopen (testOutFile.c_str(),"w");

  int it = 0;
  static int itarra[3];
  static double pq[D400][D1140][D68], pu[D400][D1140][D68], fpol[D400][D1140][D68],
    flux[D400][D1141][D68],gammin[D400][D1141],gammax[D400][D1141],
    xcell[D1141],gmax0[D400][D1141],egam[D400][D1141][44],
    phcell[D1141],rcell[D1141],ycell[D1141],cosph[D1141],
    sinph[D1141],bperp[D1141],bfield[D1141],n0[D400][D1141],
    zcell[D400][D1141],betadx[D400][D1141],betady[D400][D1141],
    betadz[D400][D1141],betad[D400][D1141],gammad[D400][D1141],
    betaux[D1141],betauy[D1141],betauz[D1141],
    betau[D1141],gammau[D1141],
    nu[D68],bx[D400][D1141],by[D400][D1141],bz[D400][D1141],
    delta[D400][D1141],enofe[D400][D1141][D44],fsync[D68],
    ididg[D1141],spsdx[BLZSIM_DIM131072],spsd[BLZSIM_DIM131072],
    gcnt[D44],igcnt[D44],fsynmd[D68][D130000],nouter[D1140],
    fsscmd[D68][D130000PADDED],fmdall[D68],deltmd[D1140],dmd[D1140],
    tlf[100],betamx[D1140],betamy[D1140],
    betamz[D1140],betamr[D1140],gamamr[D1140],bmdx[D130000],
    bmdy[D130000],bmdz[D130000],bmdtot[D130000],tlf1[D1140],
    flsync[D400][D1141][D68],flcomp[D400][D1141][D68],absorb[D68][D130000],
    fsync2[D68],cosmr[D1140],alphmd[D68][D130000],dustii[D451][D22],dustf[D451][D22],
    flec[D400][D1141][D68],flssc[D400][D1141][D68],mdd[D130000],useed[D451],
    phots[D68],phalph[D68],seedpk[D451],
    abexmd[D68][D130000],psi[D1140],sinpsi[D1140],cospsi[D1140],
    tanpsi[D1140],tauxmd[D68],phi1[D1141],phi2[D1141],
    theta1[D1141],theta2[D1141],ididb[D1141],bfrac[D1141],
    ididbd[D1141],bfracd[D1141],phi2d[D1141],thet2d[D1141],
    phi1d[D1141],thet1d[D1141],angrot[D1141],cosrot[D1141],
    bu1x[D1141],bu1y[D1141],bu1z[D1141],bu2x[D1141],bu2y[D1141],
    bu2z[D1141],cu1x[D1141],cu1y[D1141],cu1z[D1141],angrtd[D1141],
    cosrtd[D1141],bu1xd[D1141],bu1yd[D1141],bu1zd[D1141],bu2xd[D1141],
    bu2yd[D1141],bu2zd[D1141],cu1xd[D1141],cu1yd[D1141],cu1zd[D1141],
    tdelr[D451], zsmd[D130000];

  static float alfmdc[D68][D130000],syseed[D68],scseed[D68];
  int icelmx[D1141], imax[D1141], ibreak;
  double pol,pqcum,pucum,pmean,polc,ai2, pcum,tanv0,cosv0,gamup,beta,
    thlos,opang,tanop,cosop,sinop,zeta,tanz,slos,clos,eta,tanxi,xi,betacs,
    dth1,dth2,n0ave,betamd,betarl,cosmd,tanmd,bup,ustob,tlfact,delt,
    gamb,glow,gmratl,gmratm,tloss,t1,t2,tlmin, eterm1,eterm2,glim,t2max,
    sinzps,coszps,tanzps,zetap, bd2,gm1,angm,bmparx,bmpary,bmparz,
    bmprpx,bmprpy,bmprpz,dtnth1,dtnth2,dotprd,
    betupx,betupy,betupz,betatx,betaty,betatz,thetat,
    sintht,costht,phit,btprpx,btprpy,btprpz,
    btparx,btpary,btparz,sx,sy,sz,slx,sly,slz,
    anx,any,anz,bux,buy,buz,betadd,gammdd,
    bparx,bpary,bparz,bprpx,bprpy,bprpz,bpar,bprp;
  double n0mean,restnu;
  double ssabs, sstau, srcfn, fq1, fsyn1, absrb1;
  // I think we only need one dummy variable, but matching the Fortran for now
  double dum1, dum2, dum3, dum4, dum5, dum6, dum7;
  int dstart, ig;

  // This line has no analog in the Fortran code, which has tdust 1) being filled in from the
  // input file (temzinp.txt) *and* 2) being part of a common block. Since I've created two
  // separate tdust's (one in BlzSimInput and one in BlzSimCommon), need to copy the value 
  // from BlzSimInput into BlzSimCommon. This need only be done once, and then the tdust in 
  // BlzSimCommon is what gets used/modified hereafter
  common.tdust = inp.tdust;

  // Initialize the random number generator. If bTestMode is true, the BlzRand::rand() method will get
  // numbers from a text file instead of actually generating new random numbers
  initRandFromTime(bTestMode);

  // Do some execution timing. For now, we want to time how long the setup takes (pre-timeloop), and the timeloop itself.
  time_t tSetupStart, tSetupEnd, tTimeloopStart, tTimeloopEnd;
  time(&tSetupStart);
  cout << "% Setup time start: " << ctime(&tSetupStart);

  const int ICELLS = 400; //  Can lower this for testing, but Marscher currently has this at 400. #cells along the axial direction. Needs to be increased to 200 eventually
  const int NEND = inp.nend; // NEND cells along each side of the hexagon
  const int JCELLS = 3*NEND*(NEND-1)+1; // JCELLS cells along the transverse direction (perp to axial dir)
  int  ancol = 2*NEND-1;
  double rbound = ancol*inp.rsize; // inp is the BlzSimInput object passed into this method
  // zrat later defines zsize = zrat*rsize/tanz
  // If zrat is decreased, size of dust arrays needs to be increase above 451
  const double ZRAT = 0.2;
  const int NZID = NEND/ZRAT;
  double amppsd = 6.0;
  const int NDIM = BLZSIM_DIM131072;
  const int ANDIM = NDIM;
  const int MDMD = 41820;
  const int MDMAX = 120000; // Meant to be a bit smaller than D130000
  // An SED will be printed out every ispecs time steps
  const int ISPECS = 1;
  int ispec = 1;
  // Set up frequencies
  int inu;
  nu[0] = 1e10;
  nu[1] = 1.5e10;
  nu[2] = 4.3e10;
  nu[3] = 8.6e10;
  nu[4] = 1e11; // ::pow(10.0, 10.0+0.25*(4)); // 4 = 5 - 1
  nu[5] = 2.3e11;
  for(inu=7; inu<=D68; inu++) {
    phots[inu-1] = 0.0;
    nu[inu-1] = ::pow(10.0, 10.0+0.25*(inu-1));
  }
  // alpha = (s-1)/2 is the underlying spectral index, where s = slope of electron E dist.
  // gmaxmn,gmaxmx are the min/max values of initial gamma_max of
  // electrons in a cell
  double gmaxmx = inp.gmrat * inp.gmaxmn;
  double gmrat_original = inp.gmrat;
  // gmin is the minimum value of initial gamma of electrons
  // 2p is the slope of the volume vs. initial gamma-max law over all cells,
  // V = V0*gamma_max**(-2p)
  double pexp = -2.0 * inp.p;
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
  // uratio is the user-specified ratio of energy density of electrons to that of mag. field
  // set the normalization of the electron energy distribution accordingly:
  n0ave = inp.uratio*inp.bave*inp.bave*(sen-2.0)/(8.0*PI*EMC2)/(::pow(inp.gmin, 2.0-sen) - ::pow(inp.gmaxmn,2.0-sen));
  // Line of sight
  thlos = inp.thlos/DEG_PER_RAD;
  slos = sin(thlos);
  clos = cos(thlos);
  common.sinz = sin(zeta);
  common.cosz = cos(zeta);
  tanz = common.sinz/common.cosz;
  double betaup2 = inp.betaup*inp.betaup, betat2 = inp.betat*inp.betat;
  gamup = 1.0/sqrt(1.0 - betaup2);

  // Compression ratio of Mach disk. Ultra-relativistic eq. of state assumed, so compression ratio is that
  // given by Hughes, Aller, & Aller (1989, ApJ, 341, 54)
  double etac = sqrt(8.0*::pow(gamup, 4) - 17.0*gamup*gamup + 9.0)/gamup;
  // Speed downstream of conical shock if turbulent velocity is ignored;
  // for setting cell length zsize and for first estimate of time delay
  const double BETADD_MIN = 0.57735;
  //betadd = sqrt(::pow(1.0-::pow(inp.betaup*cosz,2),2)+9.0*::pow(inp.betaup*inp.betaup*cosz*sinz,2))/(3.0*inp.betaup*sinz);
  sx = -common.sinz;
  sy = 0.0;
  sz = common.cosz;
  vdcalc(0,0, inp.betaup, sx, sy, sz, &anx, &any, &anz, &betadd, &gammdd, &eta);
  if((betadd <= BETADD_MIN) || (betadd >= inp.betaup)) { // go to 7891
    // write(6,9891)betaup,betadd,sinz,cosz
    printf("Program halted: betadd > sound speed %12.5E%12.5E%12.5E%12.5E\n", inp.betaup, betadd, common.sinz, common.cosz);
    exit(0); // go to 9000;
  }

  common.betd = betadd;  // 7891
  common.gamd = gammdd;
  // Length of a cylindrical cell in pc
  double zsize = 0.2*inp.rsize/tanz;  // Can set 0.2 to 2.0 for testing
  double rsize2 = inp.rsize*inp.rsize;
  double volc = PI*inp.rsize*inp.rsize*zsize;
  // Length and volume of cell in plasma proper frame
  double zsizep= zsize/common.gamd;
  double volcp = volc/common.gamd;
  double svmd = sqrt(inp.vmd);
  double delobs = 1.0/(common.gamd*(1.0-common.betd*clos));
  dstart = MDMD+ICELLS*delobs;
  // Time step in observer's frame in days
  double dtfact = (1.0-common.betd*clos)/common.betd;
  double dtime = 1190.0*zsize*dtfact*common.zred1;
  int itlast = BlzMath::toFortranInt(ndays/dtime); // time "index" of the last timestep (quit when it is >= itlast)
  BlzLog::warnScalar("itlast", itlast);
  int mdrang = BlzMath::toFortranInt(0.5*(1.0/dtfact+1.0));
  //int ip0 = randObj.rand(0) * 5000;
  int ip0 = 100;
  // Distance of shock from axis and apex of conical jet
  tanop = tan(opang);
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
  // Area filling factor of hot dust, from IR luminosity. 1.05E8 = 1E45/(parsec**2)
  double filld = 1.05e8*inp.ldust/(5.67e-5*2.0*PI*PI*inp.dtrad*inp.dtdist*::pow(inp.tdust,4));
  const int NENDZRAT = BlzMath::toFortranInt(NEND/ZRAT);

  for(id=1; id<=ICELLS+(NENDZRAT); id++) { // do 333 id=1,(icells+nend/zrat)
    zdist = inp.zdist0+(id-NENDZRAT)*zsize+zshock;
    // Calculate min & max angles of dust torus in plasma frame
    double dphi1=asin(inp.dtrad/sqrt(zdist*zdist+inp.dtdist*inp.dtdist));
    double dphi2=atan(inp.dtdist/zdist);
    dth1 = dphi2 - dphi1;
    dth2 = dphi2 + dphi1;
    common.csth1 = cos(dth1);
    common.csth2 = cos(dth2);
    common.dcsth1 = -(common.csth1-betadd)/(1.0-common.csth1*betadd);
    common.dcsth2 = -(common.csth2-betadd)/(1.0-common.csth2*betadd);
    // Doppler factor of dust torus emission in frame of cells
    // used to estimate frequency of peak intensity in plasma frame
    double tdel = gammdd*(1.0-betadd*0.5*(common.csth1+common.csth2));
    tdelr[id-1] = tdel;
    common.dsang = dth2 - dth1;
    // Calculate seed photon field from dust emission in plasma frame
    // Peak frequency of dust thermal emission for part of torus closest to shock
    seedpk[id-1] = 5.88e10*common.tdust*tdel;
    // Use this to set the frequency array of the seed photons
    for(inu=1; inu<=D22; inu++) { // do 3 inu=1,22
      common.dustnu[inu-1] = seedpk[id-1]*::pow(10.0,-1.4+(inu-1)*0.1);
      common.dusti[inu-1] = filld * seedph(common.dustnu[inu-1]);
      if(common.dusti[inu-1] < 1e-20)
        common.dusti[inu-1] = 1e-20;
      dustf[id-1][inu-1] = common.dustnu[inu-1];
      dustii[id-1][inu-1] = common.dusti[inu-1];
    }
  } // for(id=1; id<=ICELLS+NEND; id++) 333 continue

  for(id=1; id<=ICELLS+NZID; id++) { // do 336 id=1,icells+nzid
    double dflux=0.0;
    for(inu=2; inu<=D22; inu++) { // do 6 inu=2,22
      if((dustii[id-1][inu-2]<=0.0) || (dustii[id-1][inu-1]<=0)) {
        dflux = dflux + 0.5*(dustii[id-1][inu-1]+dustii[id-1][inu-2]) * (common.dustnu[inu-1]-common.dustnu[inu-2]);
      }
      else {
        double a = log10(dustii[id-1][inu-1]/dustii[id-1][inu-2])/log10(common.dustnu[inu-1]/common.dustnu[inu-2]);
        dflux = dflux + dustii[id-1][inu-2]/(a+1.0)*common.dustnu[inu-2]*(::pow(common.dustnu[inu-1]/common.dustnu[inu-2], a+1.0) - 1.0);
      }
    } // 6 continue
    useed[id-1] = FOURPI * dflux/C_CM_PER_SEC;
  } // 336 continue
  
  gmratl = log10(gmaxmx/inp.gmin)/40.0;

  double tinc = dtime;
  psdsim(NDIM, -inp.psdslp, -inp.psdslp, 1.0, tinc, spsdx); 
  // if(nTestOut==6) {
  //   int ni;
  //   for(ni=1; ni<=NDIM; ni++) {
  //     fprintf(pfTestOut, FORMAT_ARRAY_FLOAT, ni, spsdx[ni-1]);
  //   }
  //   fclose(pfTestOut);
  //   exit(0);
  // }
  double psdsum = 0.0;
  double psdsig = 0.0;
  int ip;
  //int ipulse = dstart+100+ip0;
  double spexp=1.0/(0.5*expon+1.0);

  if(nTestOut==6) {
    fprintf(pfTestOut, FORMAT_STRING_FLOAT, "spexp", 0, spexp);
  }

  for(ip=1; ip<=NDIM; ip++) {   // do 4997 ip=1,ndim
    // Normalize spsd by standard deviation
    psdsig = psdsig + spsdx[ip-1]*spsdx[ip-1]/(ANDIM-1.0);
  } // 4997 continue

  if(nTestOut==6) {
    fprintf(pfTestOut, FORMAT_STRING_FLOAT, "psdsig^2", 0, psdsig);
  }

  psdsig = sqrt(psdsig);

  if(nTestOut==6) {
    fprintf(pfTestOut, FORMAT_STRING_FLOAT, "psdsig", 0, psdsig);
  }

  const int IPLIM=129000; // just for limiting test output

  for(ip=1; ip<=NDIM; ip++) { // do 4998 ip=1,ndim
    // Next section is only for testing purposes. If using it, comment out call psdsim() above
    // spsdx[ip-1]=1.0
    // Add pulse of high energy density
    // double expsd=3.0-(ipulse-10-ip)/4.0;
    // if(ip > (ipulse-10)) 
    //   expsd=3.0+(ipulse-10-ip)/4.0;
    // if(ip >= (ipulse-20) && ip <= (ipulse))
    //   spsdx[ip-1] = exp(0.4*expsd)*spsdx[ip-1];
    // expsd=3.0-(ipulse+5-ip)/1.5;
    // if(ip > (ipulse+5))
    //   expsd=3.0+(ipulse+5-ip)/1.5;
    // if(ip >= (ipulse) && ip <= (ipulse+10))
    //   spsdx[ip-1] = exp(1.5*expsd)*spsdx[ip-1];
    // expsd=3.0-(ipulse+60-ip)/4.0;
    // if(ip > (ipulse+60))
    //   expsd=3.0+(ipulse+60-ip)/4.0;
    // if(ip >= (ipulse+50) && ip <= (ipulse+70))
    //   spsdx[ip-1] = exp(0.4*expsd)*spsdx[ip-1];
    // expsd=3.0-(ipulse+95-ip)/4.0;
    // if(ip > ipulse+95)
    //   expsd=3.0+(ipulse+95-ip)/4.0;
    // if(ip >= (ipulse+85) && ip <= (ipulse+105))
    //   spsdx[ip-1]=exp(0.35*expsd)*spsdx[ip-1];
    
    // Normalize spsd by standard deviation and take exponential of result
    // to get amplitude of flux variation
    spsd[ip-1] = exp(amppsd*spsdx[ip-1]/psdsig);
    if((nTestOut==6) && (ip>IPLIM)) {
      fprintf(pfTestOut, FORMAT_ARRAY_FLOAT, ip, spsd[ip-1]);
    }
    // Need to scale n0 and B by a different factor to get the desired amplitude
    spsd[ip-1] = ::pow(spsd[ip-1], spexp);
    if((nTestOut==6) && (ip>IPLIM)) {
      fprintf(pfTestOut, FORMAT_ARRAY_FLOAT, ip, spsd[ip-1]);
    }
    // Average of 10 time steps to smooth variations so that discreteness
    // of columns of cells does not cause artificial spikes of flux
    if(ip > 9)
      spsd[ip-1] = 0.1*(spsd[ip-1]+spsd[ip-1-1]+spsd[ip-1-2]+
                        spsd[ip-1-3]+spsd[ip-1-4]+spsd[ip-1-5]+spsd[ip-1-6]+
                        spsd[ip-1-7]+spsd[ip-1-8]+spsd[ip-1-9]);
    psdsum = psdsum+amppsd*spsd[ip-1]/(double)ANDIM;
    if((nTestOut==6) && (ip>IPLIM)) {
      fprintf(pfTestOut, FORMAT_STRING_FLOAT, "psdsum", ip, psdsum);
      fprintf(pfTestOut, FORMAT_STRING_FLOAT, "spsd", ip, spsd[ip-1]);
    }
  } //4998 continue

  if(nTestOut==6) {
    fclose(pfTestOut);
    exit(0);
  }

  // no 4999 here like in temz.f (it is effectively all commented out)

  // Set parameters of each cell at initial time
  // There are icells columns of cells, with JCELLS cells per column
  // The central cell is a Mach disk with low velocity; its observed radiation
  // is ignored, but it is an  important source of IC seed photons
  //write(4,6667)
  //write(5,6668)
  int i = 1, j = 0, jzero = 1, nnn, jcnt;
  for(nnn=1; nnn<=NEND-1; nnn++) { // do 999 nnn=1,(nend-1)
    for(jcnt=1; jcnt<=6*nnn; jcnt++) { // do 998 jcnt=1,(6*nnn)
      j++;
      int annn = nnn;
      double angc = 60.0/(DEG_PER_RAD*annn);
      rcell[j-1] = 2.0*inp.rsize*annn;
      xcell[j-1] = rcell[j-1]*cos(angc*(j-jzero));
      ycell[j-1] = rcell[j-1]*sin(angc*(j-jzero));
      cosph[j-1] = 1.0;
      sinph[j-1] = 0.0;
      double zcol = rcell[j-1]/tanz;
      imax[j-1] = BlzMath::round<double>((2.0*zcol/zsize) + 0.01);
      if(nTestOut==2) fprintf(pfTestOut, FORMAT2, j, imax[j-1]);
      if(imax[j-1] < 2)
        imax[j-1] = 2;
      // nouter[j-1]  =  approx. no. of cells between cell of interest and observer
      nouter[j-1] = imax[j-1];
      cosph[j-1] = xcell[j-1]/rcell[j-1];
      sinph[j-1] = ycell[j-1]/rcell[j-1];
      tanpsi[j-1] = rcell[j-1]/(zsvtex-zcol);
      psi[j-1] = atan(tanpsi[j-1]);
      cospsi[j-1] = cos(psi[j-1]);
      sinpsi[j-1] = tanpsi[j-1]*cospsi[j-1];
    } // for(jcnt=1; jcnt<=6*nnn; jcnt++)

    jzero = j + 1;
  } // for(nnn=1; nnn<=NEND-1; nnn++) 

  if(nTestOut==2) {
    fclose(pfTestOut);
    exit(0);
  }

  xcell[JCELLS-1] = 0.0;
  ycell[JCELLS-1] = 0.0;
  rcell[JCELLS-1] = 0.0;
  cosph[JCELLS-1] = 1.0;
  sinph[JCELLS-1] = 0.0;

  double zrf=zshock, xrf=0.0;

  // *** Set up Mach disk emission for earlier times ***
  // Compute time delay (no. of time steps) between when plasma passes Mach disk
  // and when it passes reference point of conical shock, in plasma frame
  double zmd = zrf - zshock;
  double idelmd = zmd/(zsize/common.betd);
  // Velocity parameters of plasma in Mach disk
  betadx[0][JCELLS-1] = 0.0;
  betady[0][JCELLS-1] = 0.0;
  betadz[0][JCELLS-1] = 1.0/3.0;
  betamd = betadz[0][JCELLS-1];
  double gammd = 1.0/sqrt(1.0 - betamd*betamd);
  betad[0][JCELLS-1] = betamd;
  gammad[0][JCELLS-1] = gammd;
  double dopref = gammdd/gammd;
  double dopref2 = dopref*dopref;

  i = 1;
  int md;
  double xrand, phi;
  double gminmd = 0.15*inp.gmaxmn;

  // Initialize some parameters for selecting magnetic field vector
  for(j=1; j<=JCELLS; j++) { // do 6129 j=1,jcells
    ididb[j-1] = 0;
    phi1[j-1] = 0.0;
    bfrac[j-1] = 0.0;
  } 

  // Determine B vector of MD assuming random magnetic field orientation
  // Randomly select magnetic field direction for every 10th cell, then interpolate
  // inside loop to get direction for intermediate cells
  // First need to set up field direction for the more downstream cell
  phi2[JCELLS-1] = TWOPI*randObj.rand(0);
  double costh = 2.0*(randObj.rand(0)-0.5);
  costh = adjustTrig(costh);
  theta2[JCELLS-1] = acos(costh);

  // Loop over MDMAX cells that pass through Mach disk
  for(md=1; md<=MDMAX-1; md++) { // do 130 md=1,(MDMAX-1)

    if(ididb[JCELLS-1] != 1) { // if(ididb(JCELLS).eq.1)go to 6130
      phi1[JCELLS-1] = TWOPI * randObj.rand(0);
      costh = 2.0*(randObj.rand(0)-0.5);
      costh = adjustTrig(costh);
      theta1[JCELLS-1] = acos(costh);
      double sint1 = sin(theta1[JCELLS-1]);
      double sint2 = sin(theta2[JCELLS-1]);
      bu1x[JCELLS-1] = sint1*cos(phi1[JCELLS-1]);
      bu1y[JCELLS-1] = sint1*sin(phi1[JCELLS-1]);
      bu1z[JCELLS-1] = costh;
      bu2x[JCELLS-1] = sint2*cos(phi2[JCELLS-1]);
      bu2y[JCELLS-1] = sint2*sin(phi2[JCELLS-1]);
      bu2z[JCELLS-1] = cos(theta2[JCELLS-1]);
      cosrot[JCELLS-1] = bu1x[JCELLS-1]*bu2x[JCELLS-1]+bu1y[JCELLS-1]*bu2y[JCELLS-1]+bu1z[JCELLS-1]*bu2z[JCELLS-1];
      cosrot[JCELLS-1] = adjustTrig(cosrot[JCELLS-1]);
      angrot[JCELLS-1] = acos(cosrot[JCELLS-1]);
      double xsign = randObj.rand(0)-0.5;
      xsign = xsign/abs(xsign);
      if(xsign < 0.0)
        angrot[JCELLS-1] = angrot[JCELLS-1]-TWOPI;
      cu1x[JCELLS-1] = bu1y[JCELLS-1]*bu2z[JCELLS-1]-bu1z[JCELLS-1]*bu2y[JCELLS-1];
      cu1y[JCELLS-1] = -bu1x[JCELLS-1]*bu2z[JCELLS-1]+bu1z[JCELLS-1]*bu2x[JCELLS-1];
      cu1z[JCELLS-1] = bu1x[JCELLS-1]*bu2y[JCELLS-1]-bu1y[JCELLS-1]*bu2x[JCELLS-1];
      ididb[JCELLS-1] = 1;
    }

    BlzMath::vecRot(bu1x[JCELLS-1],bu1y[JCELLS-1],bu1z[JCELLS-1],
                    cu1x[JCELLS-1],cu1y[JCELLS-1],cu1z[JCELLS-1],  
                    (bfrac[JCELLS-1]*angrot[JCELLS-1]),
                    &bux,&buy,&buz);
 
    bfrac[JCELLS-1] = bfrac[JCELLS-1]+0.1;

    // Use .99 instead of 1 to deal with precision/error issues. Really should just
    // use an integer counter for this ("every 10th cell") but just trying to match the
    // Fortran code at this point (MSV 8/10/2012)
    if(bfrac[JCELLS-1] > .99) { // if(bfrac[jcells-1].le.1.0)go to 6131
      theta2[JCELLS-1] = theta1[JCELLS-1];
      phi2[JCELLS-1] = phi1[JCELLS-1];
      bfrac[JCELLS-1] = 0.0;
      ididb[JCELLS-1] = 0;
    } // 6131 continue

    // Compute B field components downstream of Mach disk shock
    int idelay = BlzMath::toFortranInt(dstart-MDMD+md-idelmd);
    //if((idelay<1) || (idelay>(NDIM-ip0)))
    // ,  write(6,9299)idelay,i,jcells,md,idelmd,ip0,MDMAX,dstart,
    // ,  zshock,delobs,gamd,betd,clos
    //double bavg = inp.bave*sqrt(spsd[idelay+ip0-1]);
    double bavg = inp.bave;
    n0mean = n0ave*spsd[idelay+ip0-1];
    j = JCELLS;
    n0[i-1][j-1] = etac*n0mean;
    bx[i-1][j-1] = bavg*bux*etac;
    by[i-1][j-1] = bavg*buy*etac;
    bz[i-1][j-1] = bavg*buz;
    bfield[j-1] = BlzMath::mag(bx[i-1][j-1], by[i-1][j-1], bz[i-1][j-1]);
    common.bfld = bfield[j-1];
    bmdx[md-1] = bx[i-1][j-1];
    bmdy[md-1] = by[i-1][j-1];
    bmdz[md-1] = bz[i-1][j-1];
    bmdtot[md-1] = common.bfld;
    // Calculate the initial maximum electron energy in the Mach disk
    // double ball2 = BlzMath::magSquared(bux,buy,buz); // Uncomment this if uncommenting usage below
    // Next line relates this to direction of upstream B field relative to shock
    // gmax0[i-1][j-1]=gmaxmx*buz*buz/ball2
    //if(gmax0[i-1][j-1] < inp.gmaxmn) 
   //  gmax0[i-1][[j-1]=inp.gmaxmn;
    // Next 3 lines assume a power-law distribution of gmax0, unrelated
    //     to direction of B field
    // See http://mathworld.wolfram.com/RandomNumber.html for choosing
    //  random numbers from a power-law distribution
    // xrand = randObj.rand(0);
    // if(pexp == 0.0)
    //   gmax0[i-1][j-1] = gmaxmx;
    // if(pexp < 0.0)
    //  gmax0[i-1][j-1] = ::pow((::pow(gmaxmx,pexp)-::pow(inp.gmaxmn,pexp))*xrand+::pow(inp.gmaxmn,pexp),1.0/pexp);
    gmax0[i-1][j-1] = inp.gmaxmn;
    gminmd = inp.gmin;
    // gminmd=1800.0
    // gmax0(i,j)=1.4*gminmd
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

    // Effective length of Mach disk zsmd = loss time of highest-E
    // electrons (to energy glow where radiation is negligible) times flow velocity
    glow = 10.0;
    delt = (gmax0[i-1][j-1]-glow)/(tlfact*glow*gmax0[i-1][j-1]);
    zsmd[md-1] = delt/(SEC_PER_YEAR*3.26)*gammd*betamd;

    // Use 0.99 instead of 1 to avoid singularity of N(gamma) at gmax0
    inp.gmrat = 0.99*gmax0[i-1][j-1]/gminmd; // ick, Fortran is modifying the input variable
    gmratl = log10(inp.gmrat)/32.0; // TODO: copy inp.gmrat to a different variable, maybe common.gmrat
    gmratm = log10(gminmd/glow)/11.0;

    for(ie=1; ie<=D44; ie++) { // do 124 ie=1,44
      if(ie >= 12)
        egam[i-1][j-1][ie-1] = gminmd*::pow(10.0, gmratl*(ie-12));
      if(ie < 12)
        egam[i-1][j-1][ie-1] = glow*::pow(10.0, gmratm*(ie-1));
      double powTerm1=0.0, powTerm2=0.0;
      common.ggam[ie-1] = egam[i-1][j-1][ie-1];
      t2 = (gmax0[i-1][j-1] - common.ggam[ie-1])/(tlfact*common.ggam[ie-1]*gmax0[i-1][j-1]);
      enofe[i-1][j-1][ie-1] = 0.0;
      eterm1 = common.ggam[ie-1]/gminmd;
      eterm2 =1.0 - common.ggam[ie-1]*t2*tlfact;
      if(eterm2>=0.0) { // goto 5123
        if(ie>=12) {
          powTerm1 =  (1.0 - ::pow(eterm2, sen-1.0));
          enofe[i-1][j-1][ie-1] = n0[i-1][j-1]/((sen-1.0)*tlfact*::pow(common.ggam[ie-1],sen+1.0))
            * powTerm1;
        }
        if(ie<12) {
          powTerm2 = (::pow(eterm1, sen-1.0) - ::pow(eterm2, sen-1.0));
          enofe[i-1][j-1][ie-1] = n0[i-1][j-1]/((sen-1.0)*tlfact*::pow(common.ggam[ie-1], sen+1.0))*
            powTerm2;
        }
        // Divide by delt since integral is over time
        delt = zsize*SEC_PER_YEAR*3.26/(gammad[i-1][j-1]*betad[i-1][j-1]);
        enofe[i-1][j-1][ie-1] = enofe[i-1][j-1][ie-1]/delt;
        if(enofe[i-1][j-1][ie-1]<0.0)
          enofe[i-1][j-1][ie-1] = 0.0;
      }
      common.edist[ie-1] = enofe[i-1][j-1][ie-1]; // 5123
      if((nTestOut==7) && (md>119998)) {
        fprintf(pfTestOut, FORMAT_EDIST, "5123", i, j, ie, md, powTerm1, powTerm2, n0[i-1][j-1], n0mean, etac, common.edist[ie-1]);
      }
    } // for(ie=1; ie<=D44; ie++) 124 continue
    
    if(nTestOut==7) {
      continue; // we've outputed edist, just move on
    }

    common.bperpp = common.bfld;
    common.nuhi = 1;

    for(inu=1; inu<=40; inu++) { // 125  TODO: what is the significance of inu 40??
      alfmdc[inu-1][md-1] = 10.0;
      fsscmd[inu-1][md-1] = 0.0;
      common.snu[inu-1] = nu[inu-1];
      restnu = nu[inu-1];
      // Synchrotron mean intensity for SSC calculation inside Mach disk
      ssabs = 1.02e4*(sen+2.0)*akapnu(restnu)*common.bperpp/(nu[inu-1]*nu[inu-1]);
      sstau = ssabs*CM_PER_PARSEC*inp.rsize*svmd;
      fsync[inu-1]=ajnu(restnu)*common.bperpp*inp.rsize*svmd*(EMFACT*CM_PER_PARSEC);
      //if(nTestOut==3) fprintf(pfTestOut, FORMAT3_1, md, inu, fsync[inu-1]);
      if(fsync[inu-1] > 0.0) 
        common.nuhi=inu;
      if(sstau >= 0.01) {
        srcfn = fsync[inu-1]/sstau;
        if(sstau > 5.0)
          fsync[inu-1] = srcfn;
        if((sstau>0.01) && (sstau<=5.0))
          fsync[inu-1]= srcfn*(1.0-exp(-sstau));
      } // 126
      //if(nTestOut==3) fprintf(pfTestOut, FORMAT3_1, md, inu, fsync[inu-1]);
      // Now estimate synchrotron emission seen by other cells
      // Need to add a lower frequency to get spectral index of nu(1)
      if(inu <= 1) {
        fq1 = 0.98*nu[0];
        fsyn1 = ajnu(fq1/dopref)*dopref*dopref*common.bperpp/dtfact;
        absrb1 = 1.02e4*(sen+2.0)*akapnu(fq1/dopref)*common.bperpp/::pow(fq1/dopref, 2);
      } // 127 continue
      fsynmd[inu-1][md-1] = ajnu(restnu/dopref)*dopref*dopref*common.bperpp/dtfact;
      absorb[inu-1][md-1] = 1.02e4*(sen+2.0)*akapnu(restnu/dopref)*common.bperpp/(nu[inu-1]*nu[inu-1]);
      alphmd[inu-1][md-1] = 10.0;
      abexmd[inu-1][md-1] = 1.7;
      if((fsynmd[inu-1][md-1] > 0.0) && (fsyn1 > 0.0))
        alphmd[inu-1][md-1] = -log10(fsynmd[inu-1][md-1]/fsyn1)/log10(restnu/fq1);
      if((absorb[inu-1][md-1] > 0.0) && (absrb1 > 0.0))
        abexmd[inu-1][md-1] = -log10(absorb[inu-1][md-1]/absrb1)/log10(restnu/fq1);
      fq1 = restnu;
      fsyn1 = fsynmd[inu-1][md-1];
      absrb1 = absorb[inu-1][md-1];
      common.ssseed[inu-1] = fsync[inu-1];
      //if(nTestOut==3) fprintf(pfTestOut, FORMAT3_2, md, inu, fsynmd[inu-1][md-1]);
    } // for(inu=1; inu<=40; inu++)  125

    for(inu=41; inu<=D68; inu++) { // 128
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
    double fssc1 = ssc(fq1/dopref)*dopref2/dtfact;
    int tid;
    omp_set_num_threads(NUM_THREADS);
    parallelTimer.start();

    #pragma omp parallel shared(indexTracker7), private(tid, inu, restnu)
    {
      tid = omp_get_thread_num();
      int previousInu = threadIntervals7[tid][0]; // we never use the endpoint of each interval [tid][1]
      while((inu = indexTracker7.getNextIndex(previousInu)) >= 0) { // 129
        restnu = nu[inu-1];
        double anuf = restnu/dopref;
        sscTimer.start(tid);
        fsscmd[inu-1][md-1] = ssc(anuf)*dopref2/dtfact;
        sscTimer.end(tid);
        previousInu = inu;
      }
    } // end pragma parallel

    parallelTimer.end();
    indexTracker7.reset(); // mark all indices as unused, so we can use this tracker later

    for(inu=7; inu<=D68; inu++) { // 129  why inu 7?
      restnu = nu[inu-1];
      alfmdc[inu-1][md-1] = 10.0;
      if((fsscmd[inu-1][md-1] > 0.0) && (fssc1 > 0.0))
        alfmdc[inu-1][md-1] = -log10(fsscmd[inu-1][md-1]/fssc1)/log10(restnu/fq1);
      fq1 = restnu;
      fssc1 = fsscmd[inu-1][md-1];
      if(nTestOut==3) fprintf(pfTestOut, FORMAT3_3, md, inu, fsscmd[inu-1][md-1]);
    } // 129 continue

    if((nTestOut==3) && (md > 1))
      break;

  } // for(md=1; md<=MDMAX-1; md++) 130 continue

  if(nTestOut==3) {     
    cout << "% Setup parallel time: " << parallelTimer.getTotalTime() << " sec" << endl;
    cout << "% Setup ssc() time: " << sscTimer.getTotalTime() << " sec" << endl;
    cout << "% ssc() timer detail: " << sscTimer.toString() << endl;
    fclose(pfTestOut); 
    exit(0); 
  }

  time(&tSetupEnd);
  cout << "% Setup time end: " << ctime(&tSetupEnd);
  cout << "% Setup time: " << (tSetupEnd-tSetupStart)/60. << " min" << endl;
  if(nTestOut==7) { fclose(pfTestOut); exit(0); }

  //
  //     *** End Mach disk set-up ***
  // 

  //
  // This is the top-level time loop. 
  //
  // Toward the end of the loop, it is compared to itlast and if they are equal, we break out of this loop
  time(&tTimeloopStart);
  cout << "% Timeloop time start: " << ctime(&tTimeloopStart);
  sscTimer.reset();

  while(true) {

    i = 1;
    int ncells = 0;
    it = it + 1;
    md = MDMAX;

    //
    // Start loop over all cells in first layer to set up physical parameters
    //
    double ssx, ssy, ssz, bterm;

    for(j=1; j<=JCELLS; j++) { // do 80

      for(inu=1; inu<=D68; inu++) { // do 81
        fsynmd[inu-1][md-1] = 0.0;
        fsscmd[inu-1][md-1] = 0.0;
        fmdall[inu-1] = 0.0;
      } // 81 continue

      ididg[j-1] = 0;
      zcell[i-1][j-1] = zshock - (rcell[j-1]-inp.rsize)/tanz;
      int idelay = BlzMath::toFortranInt(0.5+dstart+it-1+((zrf-zcell[i-1][j-1])+(xrf-xcell[j-1])*betadd*slos/(1.0-betadd*clos))/zsize);
      // if(idelay.le.0.or.idelay.gt.(NDIM-ip0))
      //,  write(5,9222)idelay,dstart,j,md,ip0,zshock,zrf,zcell(i,j),
      //,  xrf,xcell(j),betadd,zsize // need to put this back in once we create the file ("temzpol.txt") earlier in loop
      //double bavg = inp.bave * sqrt(spsd[idelay+ip0-1]);
      double bavg = inp.bave;
      n0mean = n0ave*spsd[idelay+ip0-1];
      phcell[j-1] = atan2(sinph[j-1],cosph[j-1]);
      // Velocity vector of laminar component of pre-shock flow
      double betupx = inp.betaup*cosph[j-1]*sinpsi[j-1];
      double betupy = inp.betaup*sinph[j-1]*sinpsi[j-1];
      double betupz = inp.betaup*cospsi[j-1];
      // Velocity vector of the turbulent component of pre-shocked plasma
      phit = TWOPI*randObj.rand(0);
      costht = 2.0*(randObj.rand(0)-0.5);
      costht = adjustTrig(costht);
      thetat = acos(costht);
      sintht = sin(thetat);
      betatx = inp.betat*cos(phit)*sintht;
      betaty = inp.betat*sin(phit)*sintht;
      betatz = inp.betat*costht;
      double dotprd = betupx*betatx+betupy*betaty+betupz*betatz;
      double btparx = dotprd*betupx/betaup2;
      double btpary = dotprd*betupy/betaup2;
      double btparz = dotprd*betupz/betaup2;
      double btprpx = betatx-btparx;
      double btprpy = betaty-btpary;
      double btprpz = betatz-btparz;
      // Velocity vector of the pre-shock plasma including turbulent component
      double gamup2 = gamup*gamup;
      betaux[j-1] = (betupx+btparx+btprpx/gamup)/(1.0+dotprd);
      betauy[j-1] = (betupy+btpary+btprpy/gamup)/(1.0+dotprd);
      betauz[j-1] = (betupz+btparz+btprpz/gamup)/(1.0+dotprd);
      betau[j-1] = BlzMath::mag(betaux[j-1], betauy[j-1], betauz[j-1]);
      gammau[j-1] = 1.0/sqrt(1.0-betau[j-1]*betau[j-1]);
      // Unit vector of shock front at current position
      sx = -common.sinz*cosph[j-1];
      sy =-common.sinz*sinph[j-1];
      sz = common.cosz;
      // Velocity vector downstream of shock + compression ratio of shock
      vdcalc(betaux[j-1],betauy[j-1],betauz[j-1],sx,sy,sz,&(betadx[i-1][j-1]),
             &(betady[i-1][j-1]),&(betadz[i-1][j-1]),&(betad[i-1][j-1]),&(gammad[i-1][j-1]),&eta);
      // Determine B vector of cell assuming random magnetic field orientation
      // 11 continue (line from the Fortran code to help me orient myself :-)  )
      // For Mach disk, continue from previous calculation of B direction
      if((j!=JCELLS) || (ididb[j-1]!=1)) { // go to 6030
        if(j != JCELLS) { // go to 6028
          if(ididb[j-1] != 1) { // go to 6030
            // Randomly select magnetic field direction for every 10th cell, then interpolate
            // inside loop to get direction for intermediate cells
            // First need to set up field direction for the more downstream cell
            phi2[j-1] = TWOPI*randObj.rand(0);
            costh = 2.0*(randObj.rand(0)-0.5);
            costh = adjustTrig(costh);
            theta2[j-1] = acos(costh);
            double sint1 = sin(theta1[j-1]);
            double sint2 = sin(theta2[j-1]);
            bu1x[j-1] = sint1*cos(phi1[j-1]);
            bu1y[j-1] = sint1*sin(phi1[j-1]);
            bu1z[j-1] = cos(theta1[j-1]);
            bu2x[j-1] = sint2*cos(phi2[j-1]);
            bu2y[j-1] = sint2*sin(phi2[j-1]);
            bu2z[j-1] = costh;
            cosrot[j-1] = bu1x[j-1]*bu2x[j-1]+bu1y[j-1]*bu2y[j-1]+bu1z[j-1]*bu2z[j-1];
            cosrot[j-1] = adjustTrig(cosrot[j-1]);
            angrot[j-1] = acos(cosrot[j-1]);
            double xsign = randObj.rand(0)-0.5;
            xsign = xsign/abs(xsign);
            if(xsign < 0.0)
              angrot[j-1] = angrot[j-1]-TWOPI;
            cu1x[j-1] = bu1y[j-1]*bu2z[j-1]-bu1z[j-1]*bu2y[j-1];
            cu1y[j-1] = -bu1x[j-1]*bu2z[j-1]+bu1z[j-1]*bu2x[j-1];
            cu1z[j-1] = bu1x[j-1]*bu2y[j-1]-bu1y[j-1]*bu2x[j-1];
            ididb[j-1] = 1;
            if(phi1[j-1] == 0) {
              // Now set up field direction for cell just crossing the shock
              phi1[j-1] = TWOPI*randObj.rand(0);
              costh = 2.0*(randObj.rand(0)-0.5);
              costh = adjustTrig(costh);
              theta1[j-1] = acos(costh);
              sint1 = sin(theta1[j-1]);
              sint2 = sin(theta2[j-1]);
              bu1x[j-1] = sint1*cos(phi1[j-1]);
              bu1y[j-1] = sint1*sin(phi1[j-1]);
              bu1z[j-1] = costh;
              bu2x[j-1] = sint2*cos(phi2[j-1]);
              bu2y[j-1] = sint2*sin(phi2[j-1]);
              bu2z[j-1] = cos(theta2[j-1]);
              cosrot[j-1] = bu1x[j-1]*bu2x[j-1]+bu1y[j-1]*bu2y[j-1]+bu1z[j-1]*bu2z[j-1];
              cosrot[j-1] = adjustTrig(cosrot[j-1]);
              angrot[j-1] = acos(cosrot[j-1]);
              xsign = randObj.rand(0)-0.5;
              xsign = xsign/abs(xsign);
              if(xsign < 0.0)
                angrot[j-1] = angrot[j-1]-TWOPI;
              cu1x[j-1] = bu1y[j-1]*bu2z[j-1]-bu1z[j-1]*bu2y[j-1];
              cu1y[j-1] = -bu1x[j-1]*bu2z[j-1]+bu1z[j-1]*bu2x[j-1];
              cu1z[j-1] = bu1x[j-1]*bu2y[j-1]-bu1y[j-1]*bu2x[j-1];
            }
          } 
        }
        else { /* j==JCELLS) */
          // Now set up field direction for cell just crossing the shock
          phi1[j-1] = TWOPI*randObj.rand(0);
          costh = 2.0*(randObj.rand(0)-0.5);
          costh = adjustTrig(costh);
          theta1[j-1] = acos(costh);
          double sint1 = sin(theta1[j-1]);
          double sint2 = sin(theta2[j-1]);
          bu1x[j-1] = sint1*cos(phi1[j-1]);
          bu1y[j-1] = sint1*sin(phi1[j-1]);
          bu1z[j-1] = costh;
          bu2x[j-1] = sint2*cos(phi2[j-1]);
          bu2y[j-1] = sint2*sin(phi2[j-1]);
          bu2z[j-1] = cos(theta2[j-1]);
          cosrot[j-1] = bu1x[j-1]*bu2x[j-1]+bu1y[j-1]*bu2y[j-1]+bu1z[j-1]*bu2z[j-1];
          cosrot[j-1] = adjustTrig(cosrot[j-1]);
          angrot[j-1] = acos(cosrot[j-1]);
          double xsign = randObj.rand(0)-0.5;
          xsign = xsign/abs(xsign);
          if(xsign < 0.0)
            angrot[j-1] = angrot[j-1]-TWOPI;
          cu1x[j-1] = bu1y[j-1]*bu2z[j-1]-bu1z[j-1]*bu2y[j-1];
          cu1y[j-1] = -bu1x[j-1]*bu2z[j-1]+bu1z[j-1]*bu2x[j-1];
          cu1z[j-1] = bu1x[j-1]*bu2y[j-1]-bu1y[j-1]*bu2x[j-1];
        }
      } //  if((j!=CELLS) || (ididb[j-1]!=1))

      // 6030
      BlzMath::vecRot(bu1x[j-1],bu1y[j-1],bu1z[j-1],cu1x[j-1],cu1y[j-1],cu1z[j-1], (bfrac[j-1]*angrot[j-1]),&bux,&buy,&buz);
      bfrac[j-1]=bfrac[j-1]+0.1;
      // Use .99 here instead of 1.0 to account for accuracy variations
      if(bfrac[j-1] > .99) {  // go to 6031  (Fortran checks bfrac against 1.0)
        theta1[j-1] = theta2[j-1];
        phi1[j-1] = phi2[j-1];
        bfrac[j-1] = 0.0;
        ididb[j-1] = 0;
      } // 6031 continue

      // Compute B field components downstream of shock in the plasma frame
      // by transforming the shock normal to the upstream plasma
      // frame, and compressing the component of B parallel to the shock
      bux = bavg*bux;
      buy = bavg*buy;
      buz = bavg*buz;

      // Unit vector of shock normal at current position
      anx = common.cosz*cosph[j-1];
      any = common.cosz*sinph[j-1];
      anz = common.sinz;

      // If Mach disk, compute separately
      if(j != JCELLS) { // go to 12
        // Calculate upstream B field components parallel + perpendicular to shock front
        bdcalc(betaux[j-1],betauy[j-1],betauz[j-1], anx, any, anz, bux, buy, buz,
               &(betau[j-1]),&(gammau[j-1]),
               &bparx,&bpary,&bparz,&bprpx,&bprpy,&bprpz,&bpar,&bprp);
        bx[i-1][j-1] = eta*bparx+bprpx;
        by[i-1][j-1] = eta*bpary+bprpy;
        bz[i-1][j-1] = eta*bparz+bprpz;
      }
      else { // if(j != JCELLS)
        // Set field of plasma in central cell, which is a Mach disk
        n0[i-1][j-1] = etac*n0mean; // 12
        bx[i-1][j-1] = bux*etac;
        by[i-1][j-1] = buy*etac;
        bz[i-1][j-1] = buz;
        bfield[j-1] = BlzMath::mag(bx[i-1][j-1], by[i-1][j-1], bz[i-1][j-1]);
        common.bfld = bfield[j-1];
        bmdx[md-1] = bx[i-1][j-1];
        bmdy[md-1] = by[i-1][j-1];
        bmdz[md-1] = bz[i-1][j-1];
        bmdtot[md-1] = common.bfld;
        // Calculate the initial maximum electron energy in the Mach disk
        // double ball2 = BlzMath::magSquared(bux,buy,buz); // uncomment this if uncomment usage below
        // Next line relates this to direction of upstream B field relative to shock
        // gmax0[i-1][j-1]=gmaxmx*buz*buz/ball2;
        // if(gmax0[i-1][j-1] < inp.gmaxmn) gmax0[i-1][j-1]=inp.gmaxmn;
        // gmax0[i-1][j-1]=inp.gmaxmn;
        // Next 3 lines assume a power-law distribution of gmax0, unrelated
        //     to direction of B field
        // See http://mathworld.wolfram.com/RandomNumber.html for choosing
        //  random numbers from a power-law distribution
        //xrand = randObj.rand(0);
        //if(pexp == 0.0)
        //  gmax0[i-1][j-1] = gmaxmx;
        //if(pexp < 0.0)
        //  gmax0[i-1][j-1] = ::pow((::pow(gmaxmx,pexp)-::pow(inp.gmaxmn,pexp))*xrand+::pow(inp.gmaxmn,pexp), 1.0/pexp);
        gmax0[i-1][j-1] = inp.gmaxmn;
        gminmd = inp.gmin;
        //gminmd=1800.0;
        //gmax0[i-1][j-1] = 1.4*gminmd;
        
        // Calculate energy distribution in the Mach disk cell
        // Compute energy density of photons in Mach disk, time delayed by 1 step
        //   Ignore photons from dust torus since beaming is small in Mach disk
        int mdi = md-1;
        if(mdi <= 0)
          mdi = 1;
        double uphmd = bmdtot[mdi-1]*bmdtot[mdi-1]/(8.0*PI)*(sqrt(1.0+4.0*inp.uratio)-1.0)/2.0;
        ustob = 8.0*PI*uphmd/(bfield[j-1]*bfield[j-1]);
        tlfact = CC2*bfield[j-1]*bfield[j-1]*(1.0+ustob);

        // Effective length of Mach disk zsmd = loss time of highest-E
        // electrons (to energy glow where radiation is negligible) times flow velocity
        glow = 10.0;
        delt = (gmax0[i-1][j-1]-glow)/(tlfact*glow*gmax0[i-1][j-1]);
        zsmd[md-1] = delt/(SEC_PER_YEAR*3.26)*gammd*betamd;
        //zsmd[md-1] = 0.2*zsize;
        // Use 0.99 instead of 1 to avoid singularity of N(gamma) at gmax0
        // [ick - shouldn't be changing value of input here. but this is what the Fortran does. TODO: Fix this]
        inp.gmrat = 0.99*gmax0[i-1][j-1]/gminmd; 
        gmratl = log10(inp.gmrat)/32.0;
        gmratm = log10(gminmd/glow)/11.0;

        for(ie=1; ie<=D44; ie++) { // 1124
          if(ie >= 12)
            egam[i-1][j-1][ie-1] = gminmd*::pow(10.0, gmratl*(ie-12));
          if(ie < 12)
            egam[i-1][j-1][ie-1] = glow*::pow(10.0, gmratm*(ie-1));
        
          common.ggam[ie-1] = egam[i-1][j-1][ie-1]; // 1123 (MSV May2013: looks like label 1123 is now an orphan, same with 123)
          t2=(gmax0[i-1][j-1]-common.ggam[ie-1])/(tlfact*common.ggam[ie-1]*gmax0[i-1][j-1]);
          enofe[i-1][j-1][ie-1] = 0.0;
          eterm1 = common.ggam[ie-1]/gminmd;
          eterm2 = 1.0-common.ggam[ie-1]*t2*tlfact;
          if(eterm2>=0.0) { // go to 5124
            if(ie>=12)
              enofe[i-1][j-1][ie-1] = n0[i-1][j-1]/((sen-1.0)*tlfact*::pow(common.ggam[ie-1],sen+1.0))
                * (1.0 - ::pow(eterm2, sen-1.0));
            if(ie<12) {
              enofe[i-1][j-1][ie-1] = n0[i-1][j-1]/((sen-1.0)*tlfact*::pow(common.ggam[ie-1], sen+1.0))*
                (::pow(eterm1, sen-1.0) - ::pow(eterm2, sen-1.0));
            }
            // Divide by delt since integral is over time
            delt = zsize*SEC_PER_YEAR*3.26/(gammad[i-1][j-1]*betad[i-1][j-1]);
            enofe[i-1][j-1][ie-1] = enofe[i-1][j-1][ie-1]/delt;
            if(enofe[i-1][j-1][ie-1] < 0.0)
              enofe[i-1][j-1][ie-1] = 0.0;          
          } // if((eterm1>=0.0) && (eterm2>=0.0)) { // go to 5124

          common.edist[ie-1] = enofe[i-1][j-1][ie-1]; // 5124
          if((nTestOut==7) && (it==1)) {
            double delt_local = (ie == 1 ? delt : 0.0);
            fprintf(pfTestOut, FORMAT_EDIST, "5124", i, j, ie, delt_local, tlfact, n0[i-1][j-1], n0mean, etac, common.edist[ie-1]);
          }
        } // for(ie=1; ie<=D44; ie++)   1124 continue
    
        common.bperpp = common.bfld;
        common.nuhi = 1;

        for(inu=1; inu<=40; inu++) { // 1125
          alfmdc[inu-1][md-1] = 10.0;
          fsscmd[inu-1][md-1] = 0.0;
          common.snu[inu-1] = nu[inu-1];
          restnu = nu[inu-1];
          // Synchrotron mean intensity for SSC calculation inside Mach disk
          ssabs = 1.02e4*(sen+2.0)*akapnu(restnu)*common.bperpp/(nu[inu-1]*nu[inu-1]);
          sstau = ssabs*CM_PER_PARSEC*inp.rsize*svmd;
          // fsync[inu-1] is the optically thin intensity
          fsync[inu-1] = ajnu(restnu)*common.bperpp*inp.rsize*svmd*(EMFACT*CM_PER_PARSEC);
          if(fsync[inu-1] > 0.0)
            common.nuhi = inu;
          if(sstau >= 0.01) {
            // Optically thick source function = emis. coef./abs. coef.
            // Note that the length of the path through the source cancels out
            srcfn = fsync[inu-1]/sstau;
            if(sstau > 5.0)
              fsync[inu-1] = srcfn;
            if((sstau>0.01) && (sstau<=5.0))
              fsync[inu-1] = srcfn*(1.0-exp(-sstau));
          } // 1126 continue

          // Now estimate synchrotron emission seen by other cells
          // Need to add a lower frequency to get spectral index of nu(1)
          if(inu <= 1) {
            fq1 = 0.98*nu[0];
            fsyn1 = ajnu(fq1/dopref)*dopref2*common.bperpp/dtfact;
            absrb1 = 1.02e4*(sen+2.0)*akapnu(fq1/dopref)*common.bperpp/::pow(fq1/dopref,2);
          } //1127 continue
        
          fsynmd[inu-1][md-1] = ajnu(restnu/dopref)*dopref2*common.bperpp/dtfact;
          absorb[inu-1][md-1] = 1.02e4*(sen+2.0)*akapnu(restnu/dopref)*common.bperpp/(nu[inu-1]*nu[inu-1]);
          alphmd[inu-1][md-1] = 10.0;
          abexmd[inu-1][md-1] = 1.7;
          if((fsynmd[inu-1][md-1]>0.0) && (fsyn1>0.0))
            alphmd[inu-1][md-1] = -log10(fsynmd[inu-1][md-1]/fsyn1)/log10(restnu/fq1);
          if((absorb[inu-1][md-1]>0.0) && (absrb1>0.0))
            abexmd[inu-1][md-1] = -log10(absorb[inu-1][md-1]/absrb1)/log10(restnu/fq1);
          //  write(5,9996)restnu,fsync[inu-1],fsynmd[inu-1][md-1],sstau,
          // ,   absorb[inu-1][md-1],alphmd[inu-1][md-1],bperpp
          fq1 = restnu;
          fsyn1 = fsynmd[inu-1][md-1];
          absrb1 = absorb[inu-1][md-1];
          common.ssseed[inu-1] = fsync[inu-1];
          if(nTestOut==4) fprintf(pfTestOut, FORMAT14001, 1125, i, j, md, inu, "fsync(inu)", fsync[inu-1]);
          if(nTestOut==4) fprintf(pfTestOut, FORMAT14001, 1125, i, j, md, inu, "fsynmd(inu,md)", fsynmd[inu-1][md-1]);
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
        double fssc1;
        sscTimer.start(0);
        fssc1 = ssc(fq1/dopref)*dopref2/dtfact;
        sscTimer.end(0);
        int tid;
        omp_set_num_threads(NUM_THREADS);
        
        #pragma omp parallel shared(indexTracker7), private(tid, inu, restnu)
        {
          tid = omp_get_thread_num();
          int previousInu = threadIntervals7[tid][0];
          while((inu = indexTracker7.getNextIndex(previousInu)) >= 0) { // 1129
            restnu = nu[inu-1];
            double anuf = restnu/dopref;
            sscTimer.start(tid);
            fsscmd[inu-1][md-1] = ssc(anuf)*dopref2/dtfact;
            sscTimer.end(tid);
            previousInu = inu;
          }
        }

        indexTracker7.reset();

        for(inu=7; inu<=D68; inu++) { // do 1129 inu=7,68
          restnu = nu[inu-1];
          alfmdc[inu-1][md-1] = 10.0;
          if((fsscmd[inu-1][md-1]>0.0) && (fssc1>0.0))
            alfmdc[inu-1][md-1] = -log10(fsscmd[inu-1][md-1]/fssc1)/log10(restnu/fq1);
          //  write(5,9996)restnu,fsscmd[inu-1][md-1],alfmdc[inu-1][md-1]
          fq1 = restnu;
          fssc1 = fsscmd[inu-1][md-1];
          if(nTestOut==3) fprintf(pfTestOut, FORMAT3_3, j, inu, fsscmd[inu-1][md-1]);
        } // 1129 continue

        // Exit the for(j=1; j<=JCELLS; j++) loop
        break;
      } // if(j != JCELLS)

      // 13 continue
      // Calculate component of magnetic field that is perpendicular to
      // the aberrated line of sight in the plasma frame
      // Line-of-sight unit vector
      slx = slos;
      sly = 0.0;
      slz = clos;
      bcalc(betadx[i-1][j-1],betady[i-1][j-1],betadz[i-1][j-1],
            slx,sly,slz,bx[i-1][j-1],by[i-1][j-1],bz[i-1][j-1],betau[j-1],gammau[j-1],
            &dum1,&dum2,&dum3,&dum4,&dum5,&dum6,&dum7,&(bperp[j-1]));
      bfield[j-1] = BlzMath::mag(bx[i-1][j-1], by[i-1][j-1], bz[i-1][j-1]);
      n0[i-1][j-1] = eta*n0mean;
      // Calculate the initial maximum electron energy in each cell
      // from the ratio of B(perp. to shock) to B_total in shock frame
      gmax0[i-1][j-1] = gmaxmx*(bprp*bprp/(bpar*bpar+bprp*bprp));
      if(gmax0[i-1][j-1]<inp.gmaxmn)
        gmax0[i-1][j-1] = inp.gmaxmn;
      // Next 3 lines assume a power-law distribution of gmax0, unrelated
      // to direction of B field
      // see http://mathworld.wolfram.com/RandomNumber.html for choosing
      // random numbers from a power-law distribution
      // xrand =randObj.rand(0);
      // if(pexp == 0.0)
      //   gmax0[i-1][j-1] = gmaxmx;
      // if(pexp<0.0)
      //   gmax0[i-1][j-1] = ::pow((::pow(gmaxmx,pexp)-::pow(inp.gmaxmn,pexp))*xrand+::pow(inp.gmaxmn,pexp), 1.0/pexp);
    
      for(ig=1; ig<=(D44-1); ig++) { // do 15 ig=1,43
        if((gmax0[i-1][j-1]<=gcnt[ig]) && (gmax0[i-1][j-1]>=gcnt[ig-1]))
          igcnt[ig-1] = igcnt[ig-1]+1;
      }      

    } // for(j=1; j<=JCELLS; j++)  End loop over all cells in FIRST LAYER   80 continue

    if(nTestOut==3) {
      fclose(pfTestOut);
      exit(0);
    }

    // Start another loop over all cells in first layer to calculate the
    // electron energy distribution and the emission
    double emisco, ecflux;

    for(j=1; j<=(JCELLS-1); j++) { // do 100 j=1,JCELLS-1
      ncells++;
      //if(nTestOut==2) fprintf(pfTestOut, "ncells %8d j %5d\n", ncells, j);
      emisco = 0.0;
      ecflux = 0.0;
      zcell[i-1][j-1] = zshock-(rcell[j-1]-inp.rsize)/tanz;
      // Determine Doppler factor relative to the observer
      common.bdx = betadx[i-1][j-1];
      common.bdy = betady[i-1][j-1];
      common.bdz = betadz[i-1][j-1];
      betacs = common.bdx*slos+common.bdz*clos;
      delta[i-1][j-1] = 1.0/(gammad[i-1][j-1]*(1.0-betacs));
      ecflux = 0.0;
      common.bperpp = bperp[j-1];
      common.bfld = bfield[j-1];
      // Determine velocity of the cell plasma relative to the MD plasma
      double betadij2 = betad[i-1][j-1]*betad[i-1][j-1];
      bmparx = betamd*common.bdx*common.bdz/betadij2;
      bmpary = betamd*common.bdy*common.bdz/betadij2;
      bmparz = betamd*common.bdz*common.bdz/betadij2;
      bmprpx = -bmparx;
      bmprpy = -bmpary;
      bmprpz = betamd-bmparz;
      betamx[j-1] = (common.bdx-bmparx-bmprpx/gammad[i-1][j-1])/(1.0-betamd*common.bdz);
      betamy[j-1] = (common.bdy-bmpary-bmprpy/gammad[i-1][j-1])/(1.0-betamd*common.bdz);
      betamz[j-1] = (common.bdz-bmparz-bmprpz/gammad[i-1][j-1])/(1.0-betamd*common.bdz);
      betamr[j-1] = BlzMath::mag(betamx[j-1], betamy[j-1], betamz[j-1]);
      gamamr[j-1] = 1.0/sqrt(1.0-(betamr[j-1]*betamr[j-1]));
      double zcl = zcell[i-1][j-1];
      double zrel = zshock-0.5*zsize-zcl;
      dmd[j-1] = BlzMath::mag(zrel,rcell[j-1]);
      // Calculate angle between MD seed photon and scattered l.o.s. photon
      common.cscat = scatcs(common.bdx,common.bdy,common.bdz,clos,0.0,slos,xcell[j-1],ycell[j-1],-zrel);
      // Determine Doppler factor of the MD emission in the cell's frame
      double dmdp = BlzMath::mag(rcell[j-1],zrel/gamamr[j-1]);
      cosmd = (betamx[j-1]*xcell[j-1]+betamy[j-1]*ycell[j-1]+betamz[j-1]*zrel/gamamr[j-1])/(betamr[j-1]*dmdp);
      tanmd = sqrt(1.0/(cosmd*cosmd)-1.0);
      angm = atan(tanmd);
      tanmd = tan(2.0*atan(tan(0.5*angm)/(gamamr[j-1]*(1.0+betamr[j-1]))));
      cosmd = 1.0/sqrt(1.0+tanmd*tanmd);
      deltmd[j-1] = 1.0/(gamamr[j-1]*(1.0-betamr[j-1]*cosmd));
      // Determine time step of Mach disk seed photons from light-travel delay
      double delcor = deltmd[j-1]/dopref;
      int mdmid = BlzMath::toFortranInt(MDMD-(((zrf-zcl)*clos+(xrf-xcell[j-1])*slos)+dmd[j-1])/(dtfact*zsize));
      int md1 = mdmid-mdrang;
      if(md1 < 1) md1 = 1;
      int md2 = mdmid+mdrang;
      // write(7,9191)md1,md2,mdmid
      // 9191 format(3i10)
      if(md2 > MDMAX) md2 = MDMAX;
      if(md1 > md2) md1 = md2;
      int amdrng = md2-md1+1.0;
      for(inu=1; inu<=D68; inu++) // do 3145 inu=1,68
        fmdall[inu-1] = 0.0;

      for(md=md1; md<=md2; md++) { // do 2147 md=md1,md2

        for(inu=1; inu<=D68; inu++) { // 3146
          syseed[inu-1] = 0.0;
          scseed[inu-1] = 0.0;
          tauxmd[inu-1] = 0.0;
          common.ssseed[inu-1] = 0.0;
        }

        // Mach disk's B field in frame of cell plasma transverse to cell's
        // line of sight to Mach disk (for synchrotron calculation)
        sx = rcell[j-1]*cosph[j-1];
        sy = rcell[j-1]*sinph[j-1];
        sz = zrel;
        double bmperp;
        bcalc(betamx[j-1],betamy[j-1],betamz[j-1],
              sx,sy,sz, bmdx[md-1], bmdy[md-1], bmdz[md-1],betamr[j-1], gamamr[j-1], 
              &dum1,&dum2,&dum3,&dum4,&dum5,&dum6,&dum7,&bmperp);
        // Apply as a correction factor to previous estimate of B_perpendicular
        double bpcorr=bmperp/bmdtot[md-1];

        // 
        // Calculate emission from various processes
        //

        // Calculate synchrotron flux from Mach disk in frame of cell plasma
        common.nuhi = 1;
        for(inu=1; inu<=40; inu++) { // do 146 inu=1,40
          if((fsynmd[inu-1][md-1]!=0.0) && (alphmd[inu-1][md-1]<=9.0)) {
            //  Synchrotron mean intensity for inverse Compton calculation
            ssabs = absorb[inu-1][md-1]*::pow(bpcorr, 0.5*(sen+2.0))*::pow(delcor, abexmd[inu-1][md-1]);
            double path = 2.0*inp.rsize*svmd/(tanmd*cosmd);
            if(path > zsmd[md-1]) 
              path = zsmd[md-1];
            double ajofnu = fsynmd[inu-1][md-1]*::pow(bpcorr, 1.0+alphmd[inu-1][md-1])*::pow(delcor, 2.0+alphmd[inu-1][md-1]);
            ajofnu = ajofnu*zsmd[md-1]*(EMFACT*CM_PER_PARSEC);
            syseed[inu-1] = ajofnu;
            srcfn = ajofnu/sstau;
            if(sstau > 5.0)
              syseed[inu-1] = srcfn;
            if((sstau>0.1) && (sstau<=5.0))
              syseed[inu-1] = srcfn*(1.0-exp(-sstau));
            syseed[inu-1] = syseed[inu-1]*PI*rsize2*inp.vmd/(dmd[j-1]*dmd[j-1]);
            // Now estimate exponential attenuation of MD seed photons by
            // synchrotron self-absorption in intervening cells
            double tauexp= 1.02e4*(sen+2.0)*akapnu(nu[inu-1])*bperp[j-1]/(nu[inu-1]*nu[inu-1])*CM_PER_PARSEC*dmd[j-1];
            tauxmd[inu-1] = tauexp;
            if(tauexp > 15.0)
              syseed[inu-1] = 0.0;
            if(tauexp <= 15.0)
              syseed[inu-1] = syseed[inu-1]/exp(tauexp);
          }

          scseed[inu-1] = 0.0; // 145
          if(inu < 16)
            fmdall[inu-1] = fmdall[inu-1]+syseed[inu-1]/(amdrng/dtfact);
        } // for(inu=1; inu<=40; inu++)  146 continue

        //
        // Calculate inverse Compton flux from Mach disk in frame of cell plasma
        //
        for(inu=16; inu<=D68; inu++) { // do 147 inu=16,68
          // Inverse Compton mean intensity from Mach disk for 2nd-order
          // inverse Compton calculation in cells
          scseed[inu-1] = fsscmd[inu-1][md-1]*::pow(delcor, 2.0+alfmdc[inu-1][md-1]);
          scseed[inu-1] = scseed[inu-1]*zsmd[md-1]*PI*rsize2*(AMJY*CM_PER_PARSEC)*inp.vmd/(dmd[j-1]*dmd[j-1]);

          if(tauxmd[inu-1] > 15.0)
            scseed[inu-1] = 0.0;
          if(tauxmd[inu-1] <= 15.0)
            scseed[inu-1]=scseed[inu-1]/exp(tauxmd[inu-1]);
          fmdall[inu-1] = fmdall[inu-1] + (syseed[inu-1]+scseed[inu-1])/(amdrng/dtfact);
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
        if((fmdall[inu-1]>SMALL_FMDALL) && (fmdall[inu-2]>SMALL_FMDALL)) { // go to 148
          aaa = log10(fmdall[inu-2]/fmdall[inu-1])/log10(common.snu[inu-1]/common.snu[inu-2]);
          usdmd = usdmd+0.5/C_CM_PER_SEC/(1.0-aaa)*fmdall[inu-2]*common.snu[inu-2]*(::pow(common.snu[inu-1]/common.snu[inu-2],1.0-aaa)-1.0);
        }
      } //  149 continue

      // 150 continue
      // Calculate energy distribution in the cell
      id = BlzMath::toFortranInt((zcell[i-1][j-1]-zshock)/zsize+NZID);
      if(id < 1) id = 1;

      zdist = inp.zdist0+(id-NENDZRAT)*zsize+zshock;
      // Calculate min & max angles of dust torus in plasma frame
      double dphi1 = asin(inp.dtrad/sqrt(zdist*zdist+inp.dtdist*inp.dtdist));
      double dphi2 = atan(inp.dtdist/zdist);
      dth1 = dphi2-dphi1;
      dth2 = dphi2+dphi1;
      common.csth1 = cos(dth1);
      common.csth2 = cos(dth2);
      common.dcsth1 = -(common.csth1-betad[i-1][j-1])/(1.0-common.csth1*betad[i-1][j-1]);
      common.dcsth2 = -(common.csth2-betad[i-1][j-1])/(1.0-common.csth2*betad[i-1][j-1]);
      common.csang = 0.5*(common.dcsth1+common.dcsth2);
      // Doppler factor of dust torus emission in frame of cells, used
      /// to estimate frequency of peak intensity in plasma frame
      double tdel = gammad[i-1][j-1]*(1.0-betad[i-1][j-1]*0.5*(common.csth1+common.csth2));
      // Correct energy density of seed photons from dust for Doppler factor
      // that include turbulent component of velocity; approximates that
      // motion is directed parallel to jet axis; also include first-order
      // dependence on scattering angle in Klein-Nishina cross-section
      double useedr = useed[id-1]*::pow(tdel/tdelr[id-1],2.0)/(1.0+betad[i-1][j-1]*common.csang);
      ustob = 8.0*PI*(useedr+usdmd)/(bfield[j-1]*bfield[j-1]);
      tlfact = CC2*bfield[j-1]*bfield[j-1]*(1.0+ustob);
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

        common.ggam[ie-1] = egam[i-1][j-1][ie-1]; // 89
        tloss = (gmax0[i-1][j-1]-common.ggam[ie-1])/(tlfact*common.ggam[ie-1]*gmax0[i-1][j-1]);
        t2 = min(tloss,delt);
        enofe[i-1][j-1][ie-1] = 0.0;
        eterm1 = common.ggam[ie-1]/inp.gmin;
        eterm2 = 1.0-common.ggam[ie-1]*t2*tlfact;
        if(eterm2>=0.0) { // go to 5089
          if(ie >= 12)
            enofe[i-1][j-1][ie-1] = n0[i-1][j-1]/((sen-1.0)*tlfact*::pow(common.ggam[ie-1], sen+1.0))*(1.0-::pow(eterm2,sen-1.0));
          if(ie< 12)
            enofe[i-1][j-1][ie-1] = n0[i-1][j-1]/((sen-1.0)*tlfact*::pow(common.ggam[ie-1], sen+1.0))*(::pow(eterm1,sen-1.0)-::pow(eterm2,sen-1.0));
          // Divide by delt since integral is over time
          enofe[i-1][j-1][ie-1] = enofe[i-1][j-1][ie-1]/delt;
          if(enofe[i-1][j-1][ie-1] < 0.0)
            enofe[i-1][j-1][ie-1] = 0.0;
        }
        common.edist[ie-1] = enofe[i-1][j-1][ie-1]; // 5089
        if((nTestOut==7) && (it==1)) {
          double delt_local = (ie == 1 ? delt : 0.0);
          fprintf(pfTestOut, FORMAT_EDIST, "5089", i, j, ie, delt_local, tlfact, n0[i-1][j-1], n0mean, etac, common.edist[ie-1]);
        }
      } // for(ie=1; ie<=D44; ie++)  90 continue

      delt = zsize*3.26*SEC_PER_YEAR/(gammad[i-1][j-1]*betad[i-1][j-1]);
      gammax[i-1][j-1] = gmax0[i-1][j-1]/(1.0+tlfact*delt*gmax0[i-1][j-1]);
      gammin[i-1][j-1] = inp.gmin/(1.0+tlfact*delt*inp.gmin);

      // Start loop over observed frequencies
      for(inu=1; inu<=D22; inu++) { // do 93 inu=1,22
        common.dustnu[inu-1] = dustf[id-1][inu-1];
        double hnukt = 4.8e-11*common.dustnu[inu-1]/common.tdust;
        double hnuktr = 4.8e-11*common.dustnu[inu-1]*(tdel/tdelr[id-1])/common.tdust;
        common.dusti[inu-1] = dustii[id-1][inu-1];
        if(hnuktr > 60.0) {
          common.dusti[inu-1] = 0.0;
          continue;
        }
        
        if(common.dusti[inu-1] > 1.0e-30)
          common.dusti[inu-1] = common.dusti[inu-1]*::pow(tdel/tdelr[id-1],3.0)*(1.0-exp(hnukt))/(1.0-exp(hnuktr));
      } // 93 continue

      common.betd = betad[i-1][j-1];  // don't know if this even matters
      common.gamd = gammad[i-1][j-1]; // polcalc() uses this
      double chipol = polcalc(bfield[j-1],bx[i-1][j-1],by[i-1][j-1],bz[i-1][j-1],clos,slos);

      int inuIthin = 0;

      for(inu=1; inu<=D68; inu++) {
        // Figure out which freq that ithin would be set to 1, since the parallelization of the 195 loop 
        // means that we are no longer processing the frequencies starting from inu=1.
        ssabs = 1.02e4*(sen+2.0)*akapnu(restnu)*bperp[j-1]/(nu[inu-1]*nu[inu-1])/delta[i-1][j-1];
        ssabs = ssabs*CM_PER_PARSEC*zsize;
        ssabs = ssabs*inp.rsize/zsize;
        if(ssabs <=(0.1/ancol)) {
          inuIthin = inu;
          break;
        }
        // if this works, we can get rid of ithin (msv Oct 2013)
      }

      int tid;

      #pragma omp parallel shared(indexTracker1), private(tid, inu)
      {
        tid = omp_get_thread_num();
        int previousInu = threadIntervals1[tid][0];
        while((inu = indexTracker1.getNextIndex(previousInu)) >=0 ) { // do 95, inu=1,68
          double emold = 0.0;
          double anumin, alnumn;
          int inumin, iinu;
          double restnu = nu[inu-1]*common.zred1/delta[i-1][j-1];
          common.snu[inu-1] = nu[inu-1];
          double specin = 0.0001;
          double emisco = 0.0;
          fsync2[inu-1] = 0.0;
          double ssabs = 0.0;
          double ecflux = 0.0;
          double sscflx = 0.0;
          double fsnoab = 0.0;
          flsync[i-1][j-1][inu-1] = 0.0;
          flcomp[i-1][j-1][inu-1] = 0.0;
          flec[i-1][j-1][inu-1] = 0.0;
          flssc[i-1][j-1][inu-1] = 0.0;
          flux[i-1][j-1][inu-1] = 0.0;

          if(inu <= D44) { // 92
            emisco = ajnu(restnu)*bperp[j-1]*delta[i-1][j-1]*delta[i-1][j-1];
            fsync2[inu-1] = emisco*zsize*(EMFACT*CM_PER_PARSEC);
            fsnoab = fsync2[inu-1];

            if(inu != previousInu+1) {
              // We should get here most of the time. Unless inu happens to be the next consecutive frequency
              // from the frequency that this thread just worked on, we don't have the correct value of emold.
              // So we'll have to calculate it, even if some other thread might have already calculated it.
              double previousRestnu = nu[inu-2]*common.zred1/delta[i-1][j-1]; // previous frequency, but same cell
              double previousEmisco = ajnu(previousRestnu)*bperp[j-1]*delta[i-1][j-1]*delta[i-1][j-1];
              emold = previousEmisco*zsize*(EMFACT*CM_PER_PARSEC);
            }

            if(inu != 1) { // go to 91
              if((emold>0.0) && (fsync2[inu-1]>0.0))
                specin = log10(emold/fsync2[inu-1])/log10(nu[inu-1]/nu[inu-2]);
            } // 91 continue

            if(inu <= inuIthin) { // go to 92
              ssabs = 1.02e4*(sen+2.0)*akapnu(restnu)*bperp[j-1]/(nu[inu-1]*nu[inu-1])/delta[i-1][j-1];
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
              srcfn = fsync2[inu-1]/ssabs;
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
        
          } // if(inu <= D44) 92 continue

          flsync[i-1][j-1][inu-1] = fsync2[inu-1]*(volc/zsize)*common.zred1/(1.0e18*AMJY*inp.dgpc*inp.dgpc)*FGEOM;
          flux[i-1][j-1][inu-1] = flsync[i-1][j-1][inu-1];
          if(specin < inp.alpha)
            specin = inp.alpha;
          double poldeg = (specin+1.0)/(specin+5.0/3.0);
          if(ssabs > 1.0) {
            poldeg = 3.0/(12.0*specin+19);
          }
          fpol[i-1][j-1][inu-1] = poldeg*flsync[i-1][j-1][inu-1];
          pq[i-1][j-1][inu-1] = fpol[i-1][j-1][inu-1]*cos(2.0*chipol);
          pu[i-1][j-1][inu-1] = fpol[i-1][j-1][inu-1]*sin(2.0*chipol);
          if(restnu >= 1.0e14) {
            //spxec = 0.0001;
            //spxssc = 0.0001;
            //common.betd = betad[0][JCELLS-1]; // don't see where these two lines matter...
            //common.gamd = gammad[0][JCELLS-1]; ..and they're bad for parallelization, to boot, so commenting out for now
            double ecdust_local = ecdust(restnu);
            ecflux = ecdust_local*3.086*volc*common.zred1*delta[i-1][j-1]*delta[i-1][j-1]/(inp.dgpc*inp.dgpc)*FGEOM;
            if((nTestOut==1) && (it==1)) fprintf(pfTestOut, FORMAT1, 92, i, j, inu, ecdust_local, delta[i-1][j-1]);
            sscTimer.start(tid);
            sscflx = ssc(restnu)*3.086*volc*common.zred1*delta[i-1][j-1]*delta[i-1][j-1]/(inp.dgpc*inp.dgpc)*FGEOM;
            sscTimer.end(tid);
            double taupp = 0.0;
            if(nu[inu-1] >= 1.0e22) { // go to 99
              // Pair production opacity calculation
              // Expression for anumin includes typical interaction angle
              anumin = (1.24e20/(nu[inu-1]*common.zred1))*1.24e20*::pow(2.0*gammad[i-1][j-1], 2);
              alnumn = log10(anumin);
              inumin = BlzMath::toFortranInt((alnumn-10.0)*4+1);
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
            if((nTestOut==1) && (it==1)) fprintf(pfTestOut, FORMAT1, 99, i, j, inu, ecflux, delta);
            /*
              if((emeold>0.0) && (ecflux>0.0))
              spxec = log10(emeold/ecflux)/log10(nu[inu-1]/nu[inu-2]);
              if((emsold>0.0) && (sscflx>0.0))
              spxssc = log10(emsold/sscflx)/log10(nu[inu-1]/nu[inu-2]);
            */
          } // if(restnu >= 1.0e14) 94 continue

          emold = fsnoab; // store the value of fsync2 for this value of inu, for use with the next value of inu (i.e. inu+1)
          //emeold = ecflux;
          //emsold = sscflx;

          previousInu = inu;
        } // for(inu=1; inu<=D68; inu++)  95 continue
      } // #pragma parallel

      indexTracker1.reset();

      if(nTestOut==8) {
        for(inu=1; inu<=D68; inu++) {
          double restnu =  nu[inu-1]*common.zred1/delta[i-1][j-1];
          fprintf(pfTestOut, FORMAT_STRING_3INT_FLOAT6_FLOAT8, "flec", i, j, inu, restnu, flec[i-1][j-1][inu-1]);
          fprintf(pfTestOut, FORMAT_STRING_3INT_FLOAT6_FLOAT8, "flssc", i, j, inu, restnu, flssc[i-1][j-1][inu-1]);
        }
      }

      // 96 continue Fortran has this line here, but there's no reference to it (!)
      icelmx[j-1] = 1;

    } // for(j=1; j<=(JCELLS-1); j++) 100 continue

    if(nTestOut==8) {
      fclose(pfTestOut);
      exit(0);
    }

    // 
    // End loop over first layer of cells
    //

    int istart, iend, idelay; 

    if(ICELLS != 1) {  // if(icells.eq.1) go to 2200 (skip entire downstream loop if ICELLS is 1
      i = i+1;
      istart = i;

      // Start loop over downstream cells
      for(j=1; j<=JCELLS-1; j++) { // do 200 j=1,(JCELLS-1)

        // iend, to be computed, is the last slice of cells with energetic electrons
        iend = 1;
        ididg[j-1] =0;
        for(i=istart; i<=imax[j-1]; i++) { // do 200 i=istart,imax[j-1]
          icelmx[j-1] = i;
          ncells = ncells+1;
          //if(nTestOut==2) fprintf(pfTestOut, "1908 ncells %8d j %5d i %5d\n", ncells, j, i);
          if(it <= 1) { // if(it.gt.1)go to 110
            //
            // *** Initial set-up of downstream cells; skip after 1st time step
            //
            zcell[i-1][j-1] = (i-1)*zsize+zshock-(rcell[j-1]-inp.rsize)/tanz;
            // Set up physical parameters of downstream cells at first time step
            // Cells farther downstream contain plasma ejected earlier
            idelay = BlzMath::toFortranInt(0.5+dstart+it-1+((zrf-zcell[i-1][j-1])+(xrf-xcell[j-1])*betadd*slos/(1.0-betadd*clos))/zsize);

            if((idelay<1) || (idelay>(NDIM-ip0))) {
              //,  write(5,9223)idelay,dstart,j,md,ip0,zshock,zrf,zcell(i,j),
              //,  xrf,xcell(j),betad(i,j),zsize
            }

            //double bavg = inp.bave*sqrt(spsd[idelay+ip0-1]);
            double bavg = inp.bave;
            n0mean = n0ave*spsd[idelay+ip0-1];
            // Velocity vector of laminar component of pre-shock flow
            double betupx = inp.betaup*cosph[j-1]*sinpsi[j-1];
            double betupy = inp.betaup*sinph[j-1]*sinpsi[j-1];
            double betupz = inp.betaup*cospsi[j-1];
            // Velocity vector of the turbulent component of pre-shocked plasma
            phit = 2.0*PI*randObj.rand(0);
            double costht = 2.0*(randObj.rand(0)-0.5);
            costht = adjustTrig(costht);
            double thetat = acos(costht);
            double sintht = sin(thetat);
            double betatx = inp.betat*cos(phit)*sintht;
            double betaty = inp.betat*sin(phit)*sintht;
            double betatz = inp.betat*costht;
            double dotprd = betupx*betatx+betupy*betaty+betupz*betatz;
            double betaup2 = inp.betaup*inp.betaup;
            double btparx = dotprd*betupx/betaup2;
            double btpary = dotprd*betupy/betaup2;
            double btparz = dotprd*betupz/betaup2;
            double btprpx =betatx-btparx;
            double btprpy = betaty-btpary;
            double btprpz = betatz-btparz;
            // Velocity vector of the pre-shock plasma including turbulent component
            double gamup2 = gamup*gamup;
            betaux[j-1] = (betupx+btparx+btprpx/gamup)/(1.0+dotprd);
            betauy[j-1] = (betupy+btpary+btprpy/gamup)/(1.0+dotprd);
            betauz[j-1] = (betupz+btparz+btprpz/gamup)/(1.0+dotprd);
            betau[j-1] = BlzMath::mag(betaux[j-1], betauy[j-1], betauz[j-1]);
            gammau[j-1] = 1.0/sqrt(1.0-betau[j-1]*betau[j-1]);
            // Unit vector of shock front at current position
            sx = -common.sinz*cosph[j-1];
            sy = -common.sinz*sinph[j-1];
            sz = -common.cosz;
            // Velocity vector downstream of shock + compression ratio of shock
            vdcalc(betaux[j-1],betauy[j-1],betauz[j-1],sx,sy,sz,
                   &(betadx[i-1][j-1]),&(betady[i-1][j-1]),&(betadz[i-1][j-1]),
                   &(betad[i-1][j-1]),&(gammad[i-1][j-1]),&eta);
            double betacs = betadx[i-1][j-1]*slos+betadz[i-1][j-1]*clos;
            delta[i-1][j-1] = 1.0/(gammad[i-1][j-1]*(1.0-betacs));
            // 103 continue <-- in the Fortran but can't find a reference to it
            // Randomly select magnetic field direction for every 10th cell, then interpolate
            // inside loop to get direction for intermediate cells
            // First need to initialize values to maintain continuity with it=1,i=1 values
            if(i <= 2) { // go to 6229
              ididbd[j-1]=ididb[j-1];
              bfracd[j-1]=bfrac[j-1];
              thet1d[j-1]=theta1[j-1];
              thet2d[j-1]=theta2[j-1];
              phi1d[j-1]=phi1[j-1];
              phi2d[j-1]=phi2[j-1];
              angrtd[j-1]=angrot[j-1];
              cosrtd[j-1]=cos(angrtd[j-1]);
              bu1xd[j-1]=bu1x[j-1];
              bu1yd[j-1]=bu1y[j-1];
              bu1zd[j-1]=bu1z[j-1];
              bu2xd[j-1]=bu2x[j-1];
              bu2yd[j-1]=bu2y[j-1];
              bu2zd[j-1]=bu2z[j-1];
              cu1xd[j-1]=cu1x[j-1];
              cu1yd[j-1]=cu1y[j-1];
              cu1zd[j-1]=cu1z[j-1];
            }

            if(ididbd[j-1] != 1) { //6229 go to 6230
              phi2d[j-1]= TWOPI*randObj.rand(0);
              costh=2.0*(randObj.rand(0)-0.5);
              costh = adjustTrig(costh);
              thet2d[j-1] = acos(costh);
              if(i > 2) { // go to 6228
                double sint1 = sin(thet1d[j-1]);
                double sint2 = sin(thet2d[j-1]);
                bu1xd[j-1] = sint1*cos(phi1d[j-1]);
                bu1yd[j-1] = sint1*sin(phi1d[j-1]);
                bu1zd[j-1] = cos(thet1d[j-1]);
                bu2xd[j-1] = sint2*cos(phi2d[j-1]);
                bu2yd[j-1] = sint2*sin(phi2d[j-1]);
                bu2zd[j-1] = costh;
                cosrtd[j-1] = bu1xd[j-1]*bu2xd[j-1]+bu1yd[j-1]*bu2yd[j-1]+bu1zd[j-1]*bu2zd[j-1];
                cosrtd[j-1] = adjustTrig(cosrtd[j-1]);
                angrtd[j-1] = acos(cosrtd[j-1]);
                double xsign = randObj.rand(0)-0.5;
                xsign = xsign/abs(xsign);
                if(xsign < 0.0)
                  angrtd[j-1] = angrtd[j-1]-TWOPI;
                cu1xd[j-1] = bu1yd[j-1]*bu2zd[j-1]-bu1zd[j-1]*bu2yd[j-1];
                cu1yd[j-1] = -bu1xd[j-1]*bu2zd[j-1]+bu1zd[j-1]*bu2xd[j-1];
                cu1zd[j-1] = bu1xd[j-1]*bu2yd[j-1]-bu1yd[j-1]*bu2xd[j-1];
                ididbd[j-1] = 1;
              }
              else {
                phi1d[j-1] = TWOPI*randObj.rand(0);  // 6228
                costh = 2.0*(randObj.rand(0)-0.5);
                costh = adjustTrig(costh);
                thet1d[j-1] = acos(costh);
                double sint1 = sin(thet1d[j-1]);
                double sint2 = sin(thet2d[j-1]);
                bu1xd[j-1] = sint1*cos(phi1d[j-1]);
                bu1yd[j-1] = sint1*sin(phi1d[j-1]);
                bu1zd[j-1] = costh;
                bu2xd[j-1] = sint2*cos(phi2d[j-1]);
                bu2yd[j-1] = sint2*sin(phi2d[j-1]);
                bu2zd[j-1] = cos(thet2d[j-1]);
                cosrtd[j-1] = bu1xd[j-1]*bu2xd[j-1]+bu1yd[j-1]*bu2yd[j-1]+bu1zd[j-1]*bu2zd[j-1];
                cosrtd[j-1] = adjustTrig(cosrtd[j-1]);
                angrtd[j-1] = acos(cosrtd[j-1]);
                double xsign = randObj.rand(0)-0.5;
                xsign = xsign/abs(xsign);
                if(xsign < 0.0)
                  angrtd[j-1] = angrtd[j-1]-TWOPI;
                cu1xd[j-1] = bu1yd[j-1]*bu2zd[j-1]-bu1zd[j-1]*bu2yd[j-1];
                cu1yd[j-1] = -bu1xd[j-1]*bu2zd[j-1]+bu1zd[j-1]*bu2xd[j-1];
                cu1zd[j-1] = bu1xd[j-1]*bu2yd[j-1]-bu1yd[j-1]*bu2xd[j-1];
                ididbd[j-1] = 1;
              }
            } //  if(ididbd[j-1] != 1) {  6230

            // 6230
            BlzMath::vecRot(bu1xd[j-1],bu1yd[j-1],bu1zd[j-1],cu1xd[j-1],cu1yd[j-1],cu1zd[j-1],(bfracd[j-1]*angrtd[j-1]),
                            &bux,&buy,&buz);
            bfracd[j-1] = bfracd[j-1]+0.1;
            // Use .99 here instead of 1.0 because of precision issues
            if(bfracd[j-1] > .99) { // go to 6231
              thet1d[j-1] = thet2d[j-1];
              phi1d[j-1] = phi2d[j-1];
              bfracd[j-1] = 0.0;
              ididbd[j-1] = 0;
            }
            // 6231 continue

            // Compute B field components downstream of shock in the plasma frame
            // by transforming the shock normal to the upstream plasma
            // frame, then compressing the component of B parallel to the shock
            bux = bavg*bux;
            buy = bavg*buy;
            buz = bavg*buz;
            bup = BlzMath::mag(bux, buy, buz);
            // Calculate upstream B field components parallel + perpendicular to shock front
            // Unit vector of shock normal at current position
            anx = common.cosz*cosph[j-1];
            any = common.cosz*sinph[j-1];
            anz = common.sinz;
            bdcalc(betaux[j-1],betauy[j-1],betauz[j-1],anx,any,anz,bux,buy,buz,
                   &(betau[j-1]),&(gammau[j-1]),&bparx,&bpary,&bparz,&bprpx,&bprpy,&bprpz,&bpar,&bprp);
            bx[i-1][j-1] = eta*bparx+bprpx;
            by[i-1][j-1] = eta*bpary+bprpy;
            bz[i-1][j-1] = eta*bparz+bprpz;
            n0[i-1][j-1] = eta*n0mean;

            // Calculate the initial maximum electron energy in each cell
            // Next line relates this to direction of B field relative to shock
            gmax0[i-1][j-1] = gmaxmx*(bprp*bprp/(bpar*bpar+bprp*bprp));
            if(gmax0[i-1][j-1] < inp.gmaxmn)
              gmax0[i-1][j-1] = inp.gmaxmn;
            // Next 3 lines assume a power-law distribution of gmax0, unrelated
            // to direction of B field
            // See http://mathworld.wolfram.com/RandomNumber.html for choosing
            // random numbers from a power-law distribution
            // xrand = randObj.rand(0);
            // if(pexp == 0.0)
            //   gmax0[i-1][j-1] = gmaxmx;
            // if(pexp < 0.0)
            //   gmax0[i-1][j-1] = ::pow((::pow(gmaxmx,pexp)-::pow(inp.gmaxmn,pexp))*xrand+::pow(inp.gmaxmn, pexp), 1.0/pexp);
            for(ig=1; ig<=(D44-1); ig++) { // do 105 ig=1,43
              if((gmax0[i-1][j-1]<=gcnt[ig]) && (gmax0[i-1][j-1]>=gcnt[ig-1]))
                igcnt[ig-1] = igcnt[ig-1]+1;
            }

          } // if(it <= 1) go to 110

          //
          // Time loop resumes here
          //

          // 110 continue

          // Calculate component of magnetic field that is perpendicular to
          // the aberrated line of sight in the plasma frame
          // Line-of-sight vector in plasma frame
          slx = slos;
          sly = 0.0;
          slz = clos;
          bcalc(betadx[i-1][j-1],betady[i-1][j-1],betadz[i-1][j-1],
                slx,sly,slz, bx[i-1][j-1],by[i-1][j-1],bz[i-1][j-1],betad[i-1][j-1],gammad[i-1][j-1],
                &dum1,&dum2,&dum3,&dum4,&dum5,&dum6,&dum7,&(bperp[j-1]));
          bfield[j-1] = BlzMath::mag(bx[i-1][j-1], by[i-1][j-1], bz[i-1][j-1]);
          emisco = 0.0;
          ecflux = 0.0;
          for(inu=1; inu<=D68; inu++) { // do 111 inu=1,68
            flcomp[i-1][j-1][inu-1] = 0.0;
            flsync[i-1][j-1][inu-1] = 0.0;
            flec[i-1][j-1][inu-1] = 0.0;
            flssc[i-1][j-1][inu-1] = 0.0;
            flux[i-1][j-1][inu-1] = 0.0;
          }

          common.bdx = betadx[i-1][j-1];
          common.bdy = betady[i-1][j-1];
          common.bdz = betadz[i-1][j-1];
          delt = (i-1)*zsize*SEC_PER_YEAR*3.26/(gammad[i-1][j-1]*betad[i-1][j-1]);
          common.bperpp = bperp[j-1];
          common.bfld = bfield[j-1];
          // Determine velocity of the cell plasma relative to the MD plasma
          double betadij2 = betad[i-1][j-1]*betad[i-1][j-1];
          bmparx = betamd*common.bdx*common.bdz/betadij2;
          bmpary = betamd*common.bdy*common.bdz/betadij2;
          bmparz = betamd*common.bdz*common.bdz/betadij2;
          bmprpx = -bmparx;
          bmprpy = -bmpary;
          bmprpz = betamd-bmparz;
          betamx[j-1] = (common.bdx-bmparx-bmprpx/gammad[i-1][j-1])/(1.0-betamd*common.bdz);
          betamy[j-1] = (common.bdy-bmpary-bmprpy/gammad[i-1][j-1])/(1.0-betamd*common.bdz);
          betamz[j-1] = (common.bdz-bmparz-bmprpz/gammad[i-1][j-1])/(1.0-betamd*common.bdz);
          betamr[j-1] = BlzMath::mag(betamx[j-1], betamy[j-1], betamz[j-1]);
          gamamr[j-1] = 1.0/sqrt(1.0-betamr[j-1]*betamr[j-1]);
          double zcl = zcell[i-1][j-1];
          double zrel = zshock-0.5*zsize-zcl;
          dmd[j-1] = BlzMath::mag(zrel,rcell[j-1]);
          // Calculate angle between MD seed photon and scattered l.o.s. photon
          common.cscat = scatcs(common.bdx,common.bdy,common.bdz,clos,0.0,  slos,xcell[j-1],ycell[j-1],-zrel);
          // Determine Doppler factor of the MD emission in the cell's frame
          double dmdp = BlzMath::mag(rcell[j-1],zrel/gamamr[j-1]);
          cosmd = (betamx[j-1]*xcell[j-1]+betamy[j-1]*ycell[j-1]+betamz[j-1]*zrel/gamamr[j-1])/(betamr[j-1]*dmdp);
          tanmd = sqrt(1.0/(cosmd*cosmd)-1.0);
          angm = atan(tanmd);
          tanmd = tan(2.0*atan(tan(0.5*angm)/(gamamr[j-1]*(1.0+betamr[j-1]))));
          cosmd = 1.0/sqrt(1.0+tanmd*tanmd);
          deltmd[j-1] = 1.0/(gamamr[j-1]*(1.0-betamr[j-1]*cosmd));
          // Determine time step of Mach disk seed photons from light-travel delay
          double delcor = deltmd[j-1]/dopref;
          int mdmid = BlzMath::toFortranInt(MDMD-(((zrf-zcl)*clos+(xrf-xcell[j-1])*slos)+dmd[j-1])/(dtfact*zsize));
          int md1 = mdmid-mdrang;
          if(md1 < 1)
            md1 = 1;
          int md2 = mdmid+mdrang;
          if(md2 > MDMAX)
            md2 = MDMAX;
          if(md1 > md2)
            md1 = md2;
          int amdrng = md2-md1+1.0;
          for(inu=1; inu<=D68; inu++) // do 4145 inu=1,68
            fmdall[inu-1] = 0.0l; // 4145

          for(md=md1; md<=md2; md++) { // do 4147 md=md1,md2

            for(inu=1; inu<=D68; inu++) { // do 4146 inu=1,68
              syseed[inu-1] = 0.0;
              scseed[inu-1] = 0.0;
              common.ssseed[inu-1] = 0.0;
              tauxmd[inu-1] = 0.0;
            }

            // Mach disk's B field in frame of cell plasma transverse to cell's
            // line of sight to Mach disk (for synchrotron calculation)
            sx=rcell[j-1]*cosph[j-1];
            sy=rcell[j-1]*sinph[j-1];
            sz=zrel;
            double bmperp;
            bcalc(betamx[j-1],betamy[j-1],betamz[j-1],
                  sx,sy,sz, bmdx[md-1],bmdy[md-1],bmdz[md-1],betamr[j-1],gamamr[j-1],
                  &dum1,&dum2,&dum3,&dum4,&dum5,&dum6,&dum7,&bmperp);
            // Apply as a correction factor to previous estimate of B_perpendicular
            double bpcorr = bmperp/bmdtot[md-1];

            //
            //  Calculate emission from various processes
            //
            // Calculate synchrotron flux from Mach disk in frame of cell plasma
            //
            common.nuhi = 1;
            for(inu=1; inu<=40; inu++) { // do 1146 inu=1,40
              if((fsynmd[inu-1][md-1]!=0.0) && (alphmd[inu-1][md-1]<=9.0)) { //  go to 1145
                // Synchrotron mean intensity for inverse Compton calculation
                ssabs = absorb[inu-1][md-1]*::pow(bpcorr, 0.5*(sen+2.0))*::pow(delcor, abexmd[inu-1][md-1]);
                double path = 2.0*inp.rsize*svmd/(tanmd*cosmd);
                if(path > zsmd[md-1]) 
                  path = zsmd[md-1];
                sstau = ssabs*CM_PER_PARSEC*path;
                double ajofnu = fsynmd[inu-1][md-1]*::pow(bpcorr, 1.0+alphmd[inu-1][md-1])*::pow(delcor, 2.0+alphmd[inu-1][md-1]);
                ajofnu = ajofnu*zsmd[md-1]*(EMFACT*CM_PER_PARSEC);
                syseed[inu-1] = ajofnu;
                srcfn = ajofnu/sstau;
                if(sstau > 5.0)
                  syseed[inu-1] = srcfn;
                if((sstau>0.1) && (sstau<=5.0))
                  syseed[inu-1] = srcfn*(1.0-exp(-sstau));
                syseed[inu-1] = syseed[inu-1]*PI*rsize2*inp.vmd/(dmd[j-1]*dmd[j-1]);
                double tauexp = 1.02e4*(sen+2.0)*akapnu(nu[inu-1])*bperp[j-1]/(nu[inu-1]*nu[inu-1])*CM_PER_PARSEC*dmd[j-1];
                tauxmd[inu-1] = tauexp;
                if(tauexp > 15.0)
                  syseed[inu-1] = 0.0;
                if(tauexp<=15.0)
                  syseed[inu-1] = syseed[inu-1]/exp(tauexp);
              } // if((fsynmd[inu-1][md-1]!=0.0) && (alphmd[inu-1][md-1]<=9.0))
 
              scseed[inu-1]=0.0; // 1145
              if(inu < 16)
                fmdall[inu-1] = fmdall[inu-1]+syseed[inu-1]/(amdrng/dtfact);
            } //  1146 continue

            //
            // Calculate inverse Compton flux from Mach disk in frame of cell plasma
            //

            for(inu=16; inu<=D68; inu++) { // do 1147 inu=16,68
              // Inverse Compton mean intensity from Mach disk for 2nd-order
              // inverse Compton calculation in cells
              scseed[inu-1] = fsscmd[inu-1][md-1]*::pow(delcor,2.0+alfmdc[inu-1][md-1]);
              scseed[inu-1] = scseed[inu-1]*zsmd[md-1]*PI*rsize2*(AMJY*CM_PER_PARSEC)*inp.vmd/(dmd[j-1]*dmd[j-1]);
              if(tauxmd[inu-1] > 15.0) 
                scseed[inu-1] = 0.0;
              if(tauxmd[inu-1] <= 15.0)
                scseed[inu-1]=scseed[inu-1]/exp(tauxmd[inu-1]);
              fmdall[inu-1] = fmdall[inu-1] + (syseed[inu-1]+scseed[inu-1])/(amdrng/dtfact);
            } // 1147 continue
          } //  for(md=md1; md<=md2; md++) 4147 continue
          
          common.nuhi = 0;
          for(inu=1; inu<=D68; inu++) { // do 4148 inu=1,68
            common.ssseed[inu-1] = fmdall[inu-1];
            if(fmdall[inu-1] > 0.0)
              common.nuhi=inu;
          }

          // Calculate seed photon energy density energy loss calculation
          double usdmd = 0.0, aaa;
          for(inu=2; inu<=common.nuhi; inu++) { // do 1149 inu=2,nuhi
            // Only calculate up to Klein-Nishina limit of gmin
            double epslon = 6.63e-27*common.snu[inu-1]/(inp.gmin*EMC2);
            if(epslon > 1.0)
              break; // go to 1150
            if((fmdall[inu-1]>SMALL_FMDALL) && (fmdall[inu-2]>SMALL_FMDALL)) { // go to 1148
              aaa = log10(fmdall[inu-2]/fmdall[inu-1])/log10(common.snu[inu-1]/common.snu[inu-2]);
              usdmd = usdmd+0.5/C_CM_PER_SEC/(1.0-aaa)*fmdall[inu-2]*common.snu[inu-2]*(::pow(common.snu[inu-1]/common.snu[inu-2], 1.0-aaa)-1.0);
            }  // 1148 continue
          } //for(inu=2; inu<=common.nuhi; inu++) 1149 continue

          // 1150 continue
          id = BlzMath::toFortranInt((zcell[i-1][j-1]-zshock)/zsize+NZID);
          if(id < 1)
            id = 1;

          // Skip flux calculation for cell if gamma_max is too low to emit at lowest frequency
          // Ratio of energy density of photons emitted by hot dust + Mach disk to
          // energy density of the magnetic field
          zdist = inp.zdist0+(id-NENDZRAT)*zsize+zshock;
          // Calculate min & max angles of dust torus in plasma frame
          double dphi1 = asin(inp.dtrad/sqrt(zdist*zdist+inp.dtdist*inp.dtdist));
          double dphi2 = atan(inp.dtdist/zdist);
          dth1 = dphi2-dphi1;
          dth2 = dphi2+dphi1;
          common.csth1 = cos(dth1);
          common.csth2 = cos(dth2);
          common.dcsth1 = -(common.csth1-betad[i-1][j-1])/(1.0-common.csth1*betad[i-1][j-1]);
          common.dcsth2 = -(common.csth2-betad[i-1][j-1])/(1.0-common.csth2*betad[i-1][j-1]);
          common.csang = 0.5*(common.dcsth1+common.dcsth2);
          // Doppler factor of dust torus emission in frame of cells, used
          /// to estimate frequency of peak intensity in plasma frame
          double tdel = gammad[i-1][j-1]*(1.0-betad[i-1][j-1]*0.5*(common.csth1+common.csth2));
          // Correct energy density of seed photons from dust for Doppler factor
          // that include turbulent component of velocity; approximates that
          // motion is directed parallel to jet axis; also include first-order
          // dependence on scattering angle in Klein-Nishina cross-section
          double useedr = useed[id-1]*::pow(tdel/tdelr[id-1],2.0)/(1.0+betad[i-1][j-1]*common.csang);
          ustob = 8.0*PI*(useedr+usdmd)/(bfield[j-1]*bfield[j-1]);
           // Calculate the maximum electron energy from gmax0 and solution to equation
          // d gamma/dt = - cc2*(b**2+8*pi*useed)*gamma**2
          tlfact = CC2*bfield[j-1]*bfield[j-1]*(1.0+ustob);
          int cellno = i;
          double tlavg = tlf1[j-1]/cellno;
          tlf[i-1] = tlfact;
          int l;
          for(l=istart; l<=i; l++) //  do 112 l=istart,i
            tlavg=tlavg+tlf[l-1]/cellno; // 112
          glim = sqrt(nu[0]*common.zred1/(2.8e6*bfield[j-1]*delta[i-1][j-1]));
          if(glim < 1.0)
            glim = 1.0;

          bool bSkipRestOfColumnCells = false;

          if(gammax[i-2][j-1] <= glim) { // go to 113
            ididg[j-1] = 1; // 113
            bSkipRestOfColumnCells = true; //go to 196          
          }
          else { // > glim
            int tlim = (gmax0[i-1][j-1]-glim)/(tlavg*glim*gmax0[i-1][j-1]);
            // skip rest of column of cells if energy losses are already severe
            if(delt < tlim) { //go to 114
              if(ididg[j-1] == 1)
                bSkipRestOfColumnCells = true;
            }
            else {
              ididg[j-1] = 1; // 113
              bSkipRestOfColumnCells = true;
            }
          }

          if(!bSkipRestOfColumnCells) {
            // Hopefully the logic above correctly reproduces the spaghetti go to logic in the Fortran
            // In any case, if we got here, do the rest of the column cells
            double emold=0.0, emeold=0.0, emsold=0.0;

            // iend is the last slice of cells with energetic electrons
            iend = i;
            // Calculate energy distribution for the cell
            id = BlzMath::toFortranInt((zcell[i-1][j-1]-zshock)/zsize+NEND);
            if(id < 1) id=1;
            delt = zsize*SEC_PER_YEAR*3.26/(gammad[i-1][j-1]*betad[i-1][j-1]);
            glow = gammin[i-2][j-1]/(1.0+tlfact*gammin[i-2][j-1]*delt);
            inp.gmrat = 0.99*gammax[i-2][j-1]/gammin[i-2][j-1];
            gmratl = log10(inp.gmrat)/32.0;
            gmratm = log10(gammin[i-2][j-1]/glow)/11.0;
            t2max = i*delt;
            gamb = gammax[i-2][j-1]/(1.0+tlfact*gammax[i-2][j-1]*delt);
            ibreak = 0;
            if((gamb<=glow) || (gamb>=gammax[i-2][j-1]))
              ibreak = 1;

            for(ie=1; ie<=D44; ie++) { // do 190 ie=1,44
              egam[i-1][j-1][ie-1] = gammin[i-2][j-1]*::pow(10.0,gmratl*(ie-12));
              if(ie < 12)
                egam[i-1][j-1][ie-1] = glow*::pow(10.0, gmratm*(ie-1));
              if(ie == 1)
                egam[i-1][j-1][ie-1] = glow+0.2*(glow*::pow(10.0,gmratm)-glow);
              if((ibreak!=1) && (egam[i-1][j-1][ie-1]>=gamb)) { // go to 189
                egam[i-1][j-1][ie-1] = gamb;
                ibreak = 1;
              }

              common.ggam[ie-1]=egam[i-1][j-1][ie-1]; // 189
              enofe[i-1][j-1][ie-1] = 0.0;
              t1 = (i-1)*delt;
              tloss = (gmax0[i-1][j-1]-common.ggam[ie-1])/(tlavg*common.ggam[ie-1]*gmax0[i-1][j-1]);
              if(tloss > t1) { // go to 188
                t2 = min(tloss,t2max);
                tlmin = t2max-(inp.gmin-common.ggam[ie-1])/(tlavg*inp.gmin*common.ggam[ie-1]);
                if(ie < 12)
                  t1 = tlmin;
                eterm1 = 1.0-common.ggam[ie-1]*t1*tlfact;
                eterm2 = 1.0-common.ggam[ie-1]*t2*tlfact;
                if((eterm1>=0.0) && (eterm2>=0.0)) { // go to 188
                  enofe[i-1][j-1][ie-1] = n0[i-1][j-1]/((sen-1.0)*tlavg*::pow(common.ggam[ie-1], sen+1.0))
                    *(::pow(eterm1,sen-1.0)-::pow(eterm2,sen-1.0));
                  // Divide by cell crossing time since integral is over time
                  enofe[i-1][j-1][ie-1] = enofe[i-1][j-1][ie-1]/delt;
                  if(enofe[i-1][j-1][ie-1] < 0.0)
                    enofe[i-1][j-1][ie-1] = 0.0;
                }
              } // if(tloss > t1) go to 188

              common.edist[ie-1] = enofe[i-1][j-1][ie-1]; // 188
              if((nTestOut==7) && (it==1)) {
                double delt_local = (ie == 1 ? delt : 0.0);
                fprintf(pfTestOut, FORMAT_EDIST, "188", i, j, ie, delt_local, tlfact, n0[i-1][j-1], n0mean, eta, common.edist[ie-1]);
              }
            } // for(ie=1; ie<=D44: ie++)  190 continue

            gammax[i-1][j-1] = gammax[i-2][j-1]/(1.0+tlfact*delt*gammax[i-2][j-1]);
            gammin[i-1][j-1] = gammin[i-2][j-1]/(1.0+tlfact*delt*gammin[i-2][j-1]);
            if(gammin[i-1][j-1] < 1.0)
              gammin[i-1][j-1] = 1.0;
            if(gammax[i-1][j-1] < 2.0)
              ididg[j-1] = 1;
            if(gammax[i-1][j-1] < 2.0)
              gammax[i-1][j-1] = 2.0;
            emold = 0.0;
            emeold = 0.0;
            emsold = 0.0;
            // calculate flux in mJy from cell
            for(inu=1; inu<=D22; inu++) { // do 193 inu=1,22
              common.dustnu[inu-1] = dustf[id-1][inu-1];
              double hnukt = 4.8e-11*common.dustnu[inu-1]/common.tdust;
              double hnuktr = 4.8e-11*common.dustnu[inu-1]*(tdel/tdelr[id-1])/common.tdust;
              common.dusti[inu-1] = dustii[id-1][inu-1];
              if(hnuktr > 60.0) {
                common.dusti[inu-1] = 0.0;
                continue;
              }

              if(common.dusti[inu-1] > 1.0e-30)
                common.dusti[inu-1] = common.dusti[inu-1]*::pow(tdel/tdelr[id-1],3.0)*(1.0-exp(hnukt))/(1.0-exp(hnuktr));
            }  // 193

            common.betd = betad[i-1][j-1];  // don't know if this even matters
            common.gamd = gammad[i-1][j-1]; // polcalc() uses this
            double chipol = polcalc(bfield[j-1],bx[i-1][j-1],by[i-1][j-1],bz[i-1][j-1],clos,slos);

            int inuIthin = 0;

            for(inu=1; inu<=D68; inu++) {
              // Figure out which freq that ithin would be set to 1, since the parallelization of the 195 loop 
              // means that we are no longer processing the frequencies starting from inu=1.
              ssabs = 1.02e4*(sen+2.0)*akapnu(restnu)*bperp[j-1]/(nu[inu-1]*nu[inu-1])/delta[i-1][j-1];
              ssabs = ssabs*CM_PER_PARSEC*zsize;
              ssabs = ssabs*inp.rsize/zsize;
              if(ssabs <=(0.1/ancol)) {
                inuIthin = inu;
                break;
              }
              // if this works, we can get rid of ithin (msv Oct 2013)
            }

            int tid;

            #pragma omp parallel shared(indexTracker1), private(tid, inu)
            {
              tid = omp_get_thread_num();
              int previousInu = threadIntervals1[tid][0];
              while((inu = indexTracker1.getNextIndex(previousInu)) >= 0) { // do 195 inu=1,68
                double emold = 0.0;
                double anumin, alnumn;
                int inumin, iinu;
                double restnu = nu[inu-1]*common.zred1/delta[i-1][j-1];
                common.snu[inu-1] = nu[inu-1];
                double specin = 0.0001;
                double emisco = 0.0;
                fsync2[inu-1] = 0.0;
                double ssabs = 0.0;
                double ecflux = 0.0;
                double sscflx = 0.0;
                double fsnoab = 0.0;
                double srcfn = 0.0;

                if(inu <= D44) { // go to 192
                  emisco = ajnu(restnu)*bperp[j-1]*delta[i-1][j-1]*delta[i-1][j-1];
                  fsync2[inu-1] = emisco*zsize*(EMFACT*CM_PER_PARSEC);
                  fsnoab = fsync2[inu-1];

                  if(inu != previousInu+1) {
                    // We should get here most of the time. Unless inu happens to be the next consecutive frequency
                    // from the frequency that this thread just worked on, we don't have the correct value of emold.
                    // So we'll have to calculate it, even if some other thread might have already calculated it.
                    double previousRestnu = nu[inu-2]*common.zred1/delta[i-1][j-1]; // previous frequency, but same cell
                    double previousEmisco = ajnu(previousRestnu)*bperp[j-1]*delta[i-1][j-1]*delta[i-1][j-1];
                    emold = previousEmisco*zsize*(EMFACT*CM_PER_PARSEC);
                  }

                  if(inu != 1) {  // go to 191
                    if((emold>0.0) && (fsync2[inu-1]>0.0))
                      specin = log10(emold/fsync2[inu-1])/log10(nu[inu-1]/nu[inu-2]);
                  } // 191 continue

                  if(inu <= inuIthin){ // go to 192
                    ssabs = 1.02e4*(sen+2.0)*akapnu(restnu)*bperp[j-1]/(nu[inu-1]*nu[inu-1])/delta[i-1][j-1];
                    ssabs = ssabs*CM_PER_PARSEC*zsize;
                    // Attenuate Mach disk emission from synchrotron absorption on way to cell
                    double ssabsm = ssabs*dmd[j-1]/zsize;
                    if(ssabsm >= 10.0)
                      common.ssseed[inu-1] = 0.0;
                    if((ssabsm<10.0) && (ssabsm>0.1))
                      common.ssseed[inu-1] = common.ssseed[inu-1]/exp(ssabsm);
                    // Return to absorption within cell
                    // Use rsize instead of zsize because of aberration
                    ssabs = ssabs*inp.rsize/zsize;
                    srcfn = fsync2[inu-1]/ssabs;
                    if(ssabs > 5.0)
                      fsync2[inu-1] = srcfn;
                    if((ssabs>0.1) && (ssabs<=5.0))
                      fsync2[inu-1] = srcfn*(1.0-exp(-ssabs));
                    // Attenuation from downstream cells along l.o.s. IF significant
                    double tauexp = (nouter[j-1] - i)*ssabs;
                    if((rcell[j-1] > (0.98*rbound)) && (xcell[j-1] <= 0.0))
                      tauexp = 0.0;
                    if(thlos == 0.0)
                      tauexp = (nouter[j-1]-1) * ssabs;
                    if(tauexp > 15.0)
                      fsync2[inu-1] = 0.0;
                    else 
                      fsync2[inu-1] = fsync2[inu-1]/exp(tauexp);
                    specin = inp.alpha;
                  } // if(inu <= inuIthin) 

                } // 192 continue

                // 192 continue
                flsync[i-1][j-1][inu-1] = fsync2[inu-1]*(volc/zsize)*common.zred1/(1.0e18*AMJY*inp.dgpc*inp.dgpc)*FGEOM;
                flux[i-1][j-1][inu-1] = flsync[i-1][j-1][inu-1];
                if(specin < inp.alpha)
                  specin = inp.alpha;
                double poldeg = (specin+1.0)/(specin+5.0/3.0);
                if(ssabs > 1.0) {
                  poldeg = 3.0/(12.0*specin+19);
                }
                fpol[i-1][j-1][inu-1]=poldeg*flsync[i-1][j-1][inu-1];
                pq[i-1][j-1][inu-1]=fpol[i-1][j-1][inu-1]*cos(2.0*chipol);
                pu[i-1][j-1][inu-1]=fpol[i-1][j-1][inu-1]*sin(2.0*chipol);
                if(restnu >= 1.0e14) {
                  //spxec = 0.0001;
                  //spxssc = 0.0001;
                  //common.betd = betad[0][JCELLS-1];
                  //common.gamd = gammad[0][JCELLS-1];
                  double ecdust_local = ecdust(restnu);
                  ecflux = ecdust_local*3.086*volc*common.zred1*delta[i-1][j-1]*delta[i-1][j-1]/(inp.dgpc*inp.dgpc)*FGEOM;
                  if((nTestOut==1) && (it==1)) fprintf(pfTestOut, FORMAT1, 192, i, j, inu, ecdust_local, delta[i-1][j-1]);
                  sscTimer.start(tid);
                  sscflx = ssc(restnu)*3.086*volc*common.zred1*delta[i-1][j-1]*delta[i-1][j-1]/(inp.dgpc*inp.dgpc)*FGEOM;
                  sscTimer.end(tid);
                  double taupp = 0.0;
                  if(nu[inu-1] >= 1.0e22) { // go to 199
                    // Pair production opacity calculation
                    // Expression for anumin includes typical interaction angle
                    anumin = (1.24e20/(nu[inu-1]*common.zred1))*1.24e20*::pow(2.0*gammad[i-1][j-1], 2);
                    alnumn = log10(anumin);
                    inumin = BlzMath::toFortranInt((alnumn-10.0)*4+1);
                    if(inumin <= 40) {  // go to 199
                      if(anumin > nu[inumin-1])
                        anumin = nu[inumin-1];

                      int iinu;
                      for(iinu=inumin; iinu<=40; iinu++) { // do 197 iinu=inumin,40;
                        double bp = sqrt(1.0-anumin/nu[iinu-1]);
                        if(bp <= 0.0)
                          bp=0.001;
                        double xsecp1 = 1.25e-25*(1.0-bp*bp)*((3.0-::pow(bp,4))*log((1.0+bp)/(1.0-bp))-2.0*bp*(2.0-bp*bp));
                        bp = sqrt(1.0-anumin/nu[iinu]);
                        if(bp <= 0.0)
                          bp=0.001;
                        double xsecp2 = 1.25e-25*(1.0-bp*bp)*((3.0-::pow(bp,4))*log((1.0+bp)/(1.0-bp))-2.0*bp*(2.0-bp*bp));
                        taupp = taupp+0.5*(phots[iinu-1]*xsecp1+phots[iinu]*xsecp2)*(nu[iinu]-nu[iinu-1])*nouter[j-1]*zsize*CM_PER_PARSEC;
                      } // 197 continue

                      if(taupp >= 0.1) { // go to 199
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

                  // 199
                  flec[i-1][j-1][inu-1] = ecflux;
                  flssc[i-1][j-1][inu-1] = sscflx;
                  flcomp[i-1][j-1][inu-1] = ecflux+sscflx;
                  flux[i-1][j-1][inu-1] = flsync[i-1][j-1][inu-1]+flcomp[i-1][j-1][inu-1];
                  if((nTestOut==1) && (it==1)) fprintf(pfTestOut, FORMAT1, 199, i, j, inu, ecflux, delta[i-1][j-1]);
                  // if((emeold>0.0) && (ecflux>0.0))
                  //   spxec = log10(emeold/ecflux)/log10(nu[inu-1]/nu[inu-2]);
                  // if((emsold>0.0) && (sscflx>0.0))
                  //   spxssc = log10(emsold/sscflx)/log10(nu[inu-1]/nu[inu-2]);
                } // if(restnu >= 1.0e14) 194 continue

                emold = fsnoab; // store the value of fsync2 for this value of inu, for use with the next value of inu (i.e. inu+1)
                //emeold = ecflux;
                //emsold = sscflx;
                previousInu = inu;
              } // for(inu=1; inu<=D68; inu++) 195 continue
            } // #pragma parallel

            indexTracker1.reset();
          } // if(bSkipRestOfColumnCells)  196 continue

        } // for(i-istart; i<=imax[j-1]; i++) 200 continue

      } // for(j=1; j<=JCELLS-1; j++) 200 continue

    } // if(icells.eq.1)go to 2200   2200 continue

    // End loop over downstream cells

    if(!bOutputFilesCreated) {
      pfSpec = fopen (specFile.c_str(),"w");
      pfLc = fopen (lcFile.c_str(),"w");
      pfPol = fopen (polFile.c_str(),"w");
      bOutputFilesCreated = true;
      std::string hformat("#No. of cells on each side of hexagonal grid: %3d\n#redshift: %5.3f  Distance in Gpc: %5.3f\n#spectral index: %4.2f  filling factor exponent: %4.2f\n#mean unshocked magnetic field: %5.3f  ratio of electron to mag. energy: %8.2E\n#cell radius (pc): %6.3f\n#Min. value of gamma_max: %8.1f  ratio of max. to min. values of gamma_max: %5.1f\n#gamma_min: %6.1f\n#upstream laminar velocity: %9.5fc  upstream turbulent velocity: %9.5fc\n#shock angle: %6.3f  viewing angle: %6.3f  opening angle: %6.3f\n#Dust temperature: %6.0f  dust luminosity %5.2fx10**45 erg/s\ndistance of center of dust torus from black hole: %3.1f pc\n#radius of torus: %3.1f pc   Distance of shock from central engine: %5.2f pc\n#Energy density of seed photons in plasma frame: %9.2E\n#-Slope of PSD: %5.1f     Area of Mach disk relative to other cell%9.2E\n# Area filling factor of dust emission: %9.2E\n"); // 6665
      int nendInt = static_cast<int>(inp.nend);
      fprintf(pfSpec, hformat.c_str(), nendInt, inp.zred, inp.dgpc, inp.alpha, inp.p, inp.bave, inp.uratio, inp.rsize, inp.gmaxmn,
             gmrat_original, inp.gmin, inp.betaup, inp.betat, (inp.zeta), (thlos*DEG_PER_RAD), (inp.opang),
              inp.tdust, inp.ldust, inp.dtdist,inp.dtrad,inp.zdist0,useed[4],inp.psdslp,inp.vmd,filld);
      fprintf(pfLc, hformat.c_str(), nendInt, inp.zred, inp.dgpc, inp.alpha, inp.p, inp.bave, inp.uratio, inp.rsize, inp.gmaxmn,
             gmrat_original, inp.gmin, inp.betaup, inp.betat, (inp.zeta), (thlos*DEG_PER_RAD), (inp.opang),
              inp.tdust, inp.ldust, inp.dtdist,inp.dtrad,inp.zdist0,useed[4],inp.psdslp,inp.vmd,filld);
      fprintf(pfPol, hformat.c_str(), nendInt, inp.zred, inp.dgpc, inp.alpha, inp.p, inp.bave, inp.uratio, inp.rsize, inp.gmaxmn,
             gmrat_original, inp.gmin, inp.betaup, inp.betat, (inp.zeta), (thlos*DEG_PER_RAD), (inp.opang),
              inp.tdust, inp.ldust, inp.dtdist,inp.dtrad,inp.zdist0,useed[4],inp.psdslp,inp.vmd,filld);
      fprintf(pfLc,"#\n#   time(d)     freq(Hz)  F(Jy Hz)  alpha     F(mJy)    freq(Hz)  F(Jy Hz)  alpha     F(muJy)   freq(Hz)  F(Jy Hz)   alpha     F(nJy)    F(EC)     F(SSC)  no. of live cells\n");
      fprintf(pfPol,"#\n#     time(d)  freq(Hz)   p(%%)   chi(deg)   freq      p      chi     freq       p      chi     freq      p       chi     freq       p      chi      freq       p      chi    freq       p       chi\n");
    }

    // 
    // Write out data and check if we should do another time step
    // If so, set up for next time step and jump back up to top of the time loop
    //
    double tflold, tflux, tsflux, tcflux, tecfl, tsscfl, tsscf2, alph, qcum, ucum, pang;
    double timeo, tfsync, tfsyno, tfl19, tfl20, tfl21, tfl32, tfl33, tfl34, tfl52, tfl53, tfl54;
    double tfcomp, tfcomx, tfec, tfssc, tfl6, tfl11, pdeg8, poldeg, pang8, pdeg12, pang12, pdeg16, pang16, pdeg20, pang20;
    double pdeg24, pang24, pdeg28,pang28, pdeg32, pang32;
    double pdeg3, pang3, pdeg6,pang6;
    double alp20a, alp20b, alp33a, alp33b, alp53a, alp53b;
    double alph20, alph33, alph53;
    int iwp;

    tflold=0.0; // labeled 299, but Fortran temz.f has no reference to 299
    std::string mapHeader("#    i     j    x(mas)     y(mas)   flux(Jy)     Q(Jy)       U(Jy)      P(%)    EVPA(deg)  tot flux    qcum      ucum     betad      it\n");

    for(inu=1; inu<=D68; inu++) { // do 500 inu=1,68
      tflux = 0.0;
      tsflux = 0.0;
      tcflux = 0.0;
      tecfl = 0.0;
      tsscfl = 0.0;
      alph = 0.0;
      qcum = 0.0;
      ucum = 0.0;
      iwp = getIwp(it);
      if((iwp==1) && (inu==3)) {
        char filnam[16];
        sprintf(filnam, "%8s%4.4d%4s", "ctemzmap", it, ".txt");
        mapFile.assign("maps/");
        mapFile.append(filnam);
        pfMap = fopen(mapFile.c_str(), "w");
        fprintf(pfMap, mapHeader.c_str());
      }
      
      for(j=1; j<=(JCELLS-1); j++) { //  do 300 j=1,(JCELLS-1)
        for(i=1; i<=icelmx[j-1]; i++) { // do 300 i=1,icelmx(j)
          qcum=qcum+pq[i-1][j-1][inu-1];
          ucum=ucum+pu[i-1][j-1][inu-1];
          tflux=tflux+flux[i-1][j-1][inu-1];
          tsflux=tsflux+flsync[i-1][j-1][inu-1];
          if(nTestOut==5) fprintf(pfTestOut, FORMAT14005, it, i, j, 0, inu, "flsync", flsync[i-1][j-1][inu-1],"tsflux", tsflux);
          tcflux=tcflux+flcomp[i-1][j-1][inu-1];
          tecfl=tecfl+flec[i-1][j-1][inu-1];
          tsscfl=tsscfl+flssc[i-1][j-1][inu-1];
          if(inu == 3) {
            // Position of cell on sky, in milliarcseconds (for z=0.859)
            double xobs=(zcell[i-1][j-1]*slos+xcell[j-1]*clos)/7.7;
            double yobs=ycell[j-1]/7.7;
            double poldgg = 0.0;
            double pang = 0.0;
            if(flsync[i-1][j-1][inu-1] > 1.0e-6) { // go to 298
              if(iwp == 1) {
                double poldg = 100.0*fpol[i-1][j-1][inu-1]/flsync[i-1][j-1][inu-1];
                double pangg = 0.5*atan2(pu[i-1][j-1][inu-1],pq[i-1][j-1][inu-1])*DEG_PER_RAD;
                std::string mapFormat1("%5d%5d %11.3E%11.3E%11.3E%11.3E%11.3E%11.3E%11.3E%11.3E%11.3E%11.3E%11.3E%5d\n");
                fprintf(pfMap, mapFormat1.c_str(), i,j,xobs,yobs,(0.001*flsync[i-1][j-1][inu-1]),
                        (0.001*pq[i-1][j-1][inu-1]),(0.001*pu[i-1][j-1][inu-1]),poldg,pangg,
                        tflux,qcum,ucum,betad[i-1][j-1],it);
              } // if(iwp == 1)
            }
          } // if(inu == 3)

          // 298 continue

        } // (for i=1; i<=icelmx[j-1]; i++)
      }

      // 300 continue
      phots[inu-1]=0.0;
      phalph[inu-1]=alph+1.0;
      if((tflold!=0.0) && (tflux!=0.0)) { // go to 390
        if(inu > 1)
          alph = log10(tflold/tflux)/log10(nu[inu-1]/nu[inu-2]);
        poldeg=0.0;
        if(tflux > 0.0)
          poldeg=sqrt(qcum*qcum+ucum*ucum)/tflux;
        pang=0.5*atan2(ucum,qcum)*DEG_PER_RAD;
        tflux=0.001*nu[inu-1]*tflux;
        tsflux=0.001*nu[inu-1]*tsflux;
        tcflux=0.001*nu[inu-1]*tcflux;
        tecfl=0.001*nu[inu-1]*tecfl;
        tsscfl=0.001*nu[inu-1]*tsscfl;
        tsscf2=0.001*nu[inu-1]*tsscf2;
      }
      // Print out SED for every ispecs-th time step
      // 390 continue
      timeo=it*dtime;
      if(inu == 1) {
        //Fortran: write(3,9990) timeo;
        fprintf(pfSpec, "\n\n\nTime = %8.2f days\n  freq(Hz))  Ftot(Jy Hz)  sp. index    Fsynch         F(EC)    F(SSC-MD)\n", timeo);
      }
      // Write SED to file
      tfsync=1.0e3*tsflux/nu[inu-1];
      if(inu == 20)
        tfsyno=1.0e3*tsflux/nu[inu-1];
      //write(3,9996)nu[inu-1],tflux,alph,tfsync,tecfl,tsscfl   9996 format(1p10e12.4)
      fprintf(pfSpec,"%12.4E%12.4E%12.4E%12.4E%12.4E%12.4E\n",nu[inu-1],tflux,alph,tfsync,tecfl,tsscfl);
      if(inu == 19)tfl19=tflux; // 391, but no reference ot it 
      if(inu == 20)tfl20=tflux;
      if(inu == 21)tfl21=tflux;
      if(inu == 32)tfl32=tflux;
      if(inu == 33)tfl33=tflux;
      if(inu == 34)tfl34=tflux;
      if(inu == 33)tfcomx=1.0e6*tcflux/nu[inu-1];
      if(inu == 33)tfec=tecfl;
      if(inu == 33)tfssc=tsscfl;
      if(inu == 52)tfl52=tflux;
      if(inu == 53)tfl53=tflux;
      if(inu == 54)tfl54=tflux;
      if(inu == 53)tfcomp=1.0e9*tcflux/nu[inu-1];
      if(inu == 53)tfec=tecfl;
      if(inu == 53)tfssc=tsscfl;
      if(inu == 11)tfl11=tflux;
      if(inu == 6)tfl6=tflux;
      tflold=1000.0*tflux/nu[inu-1];
      if(inu == 3)pdeg3=100.0*poldeg;
      if(inu == 3)pang3=pang;
      if(inu == 6)pdeg6=100.0*poldeg;
      if(inu == 6)pang6=pang;
      if(inu == 8)pdeg8=100.0*poldeg;
      if(inu == 8)pang8=pang;
      if(inu == 12)pdeg12=100.0*poldeg;
      if(inu == 12)pang12=pang;
      if(inu == 16)pdeg16=100.0*poldeg;
      if(inu == 16)pang16=pang;
      if(inu == 20)pdeg20=100.0*poldeg;
      if(inu == 20)pang20=pang;
      if(inu == 24)pdeg24=100.0*poldeg;
      if(inu == 24)pang24=pang;
      if(inu == 28)pdeg28=100.0*poldeg;
      if(inu == 28)pang28=pang;
      if(inu == 32)pdeg32=100.0*poldeg;
      if(inu == 32)pang32=pang;
    } // End loop over nu,  500
    
    // 500 continue

    // If we created a map file, close it
    if((iwp==1) && (inu==3))
      fclose(pfMap);
    alp20a=0.0;
    alp20b=0.0;
    alp33a=0.0;
    alp33b=0.0;
    alp53a=0.0;
    alp53b=0.0;
    if(tfl19>0.0 && tfl20>0.0) alp20a=log10(tfl19/tfl20)/log10(nu[20-1]/nu[19-1])+1.0;
    if(tfl20>0.0 && tfl21>0.0) alp20b=log10(tfl20/tfl21)/log10(nu[21-1]/nu[20-1])+1.0;
    if(tfl32>0.0 && tfl33>0.0) alp33a=log10(tfl32/tfl33)/log10(nu[33-1]/nu[32-1])+1.0;
    if(tfl33>0.0 && tfl34>0.0) alp33b=log10(tfl33/tfl34)/log10(nu[34-1]/nu[33-1])+1.0;
    if(tfl52>0.0 && tfl53>0.0) alp53a=log10(tfl52/tfl53)/log10(nu[53-1]/nu[52-1])+1.0;
    if(tfl53>0.0 && tfl54>0.0)alp53b =log10(tfl53/tfl54)/log10(nu[54-1]/nu[53-1])+1.0;
    alph20=0.5*(alp20a+alp20b);
    alph33=0.5*(alp33a+alp33b);
    alph53=0.5*(alp53a+alp53b);
    timeo=it*dtime;
    //  Write light curve points to file
    fprintf(pfLc, "%5d%8.2f  %10.2E%10.2E%10.2E%10.2E%10.2E%10.2E%10.2E%10.2E%10.2E%10.2E%10.2E%10.2E%10.2E%10.2E%10.2E%10.2E %5d\n", 
            it,timeo,nu[20-1],tfl20,alph20,tfsyno,nu[33-1],
            tfl33,alph33,tfcomx,nu[53-1],tfl53,alph53,tfcomp,tfec,
            tfssc,nu[6-1],tfl6,ncells);
    // Write selected polarization data to file
    // 9988 format(i5,f8.2,2x,7(1pe8.2,1x,0pf7.3,1x,f8.3,1x))
    fprintf(pfPol, "%5d%8.2f  %8.2E %7.3f %8.3f %8.2E %7.3f %8.3f %8.2E %7.3f %8.3f %8.2E %7.3f %8.3f %8.2E %7.3f %8.3f %8.2E %7.3f %8.3f %8.2E %7.3f %8.3f\n",
            it,timeo,nu[3-1],pdeg3,pang3,nu[6-1],pdeg6,pang6,nu[8-1],pdeg8,pang8,nu[12-1],pdeg12,pang12,
            nu[16-1],pdeg16,pang16,nu[20-1],pdeg20,pang20,nu[24-1],pdeg24,pang24);
    
    BlzLog::warnScalar("Done with time loop ", it);

    if(it >= itlast)
      break; // end the time loop

    // Set up next time step by shifting physical conditions 1 slice down jet
    for(j=1; j<=JCELLS; j++) { //   do 598 j=1,jcells
      for(i=icelmx[j-1]+1; i>=istart; i--) { // do 598 i=icelmx(j)+1,istart,-1
        gmax0[i-1][j-1]=gmax0[i-2][j-1];
        bx[i-1][j-1]=bx[i-2][j-1];
        by[i-1][j-1]=by[i-2][j-1];
        bz[i-1][j-1]=bz[i-2][j-1];
        n0[i-1][j-1]=n0[i-2][j-1];
        betadx[i-1][j-1]=betadx[i-2][j-1];
        betady[i-1][j-1]=betady[i-2][j-1];
        betadz[i-1][j-1]=betadz[i-2][j-1];
        betad[i-1][j-1]=betad[i-2][j-1];
        gammad[i-1][j-1]=gammad[i-2][j-1];
        delta[i-1][j-1]=delta[i-2][j-1];
      }
    } // 598 continue
    ispec=ispec-1;
    if(ispec == 0)
      ispec=ISPECS+1;
    // Move cells in time array to make room for next time step
    for(inu=1; inu<=D68; inu++) { // do 599 inu=1,68
      for(md=1; md<=MDMAX-1; md++) { // do 599 md=1,(MDMAX-1)
        alphmd[inu-1][md-1]=alphmd[inu-1][md];
        alfmdc[inu-1][md-1]=alfmdc[inu-1][md];
        fsynmd[inu-1][md-1]=fsynmd[inu-1][md];
        fsscmd[inu-1][md-1]=fsscmd[inu-1][md];
        absorb[inu-1][md-1]=absorb[inu-1][md];
      }
    } //continue
     
    // go to 9
    // Loop should be terminated by the if(it==itlast) check above
  } // while(true) Top-level time loop

  time(&tTimeloopEnd);
  cout << "% Timeloop time end: " << ctime(&tTimeloopEnd);
  cout << "% Timeloop time: " << (tTimeloopEnd-tTimeloopStart)/60. << " min" << endl;
  cout << "% Timeloop ssc() time: " << sscTimer.getTotalTime()/60. << " min" << endl;
  BlzLog::warn("Done with top level time loop");
  BlzLog::warn("Closing output files...");
  fclose(pfSpec);
  fclose(pfLc);
  fclose(pfPol);
  if(nTestOut>0) fclose(pfTestOut);
  BlzLog::warn("Done closing output files");

} // end run() method
