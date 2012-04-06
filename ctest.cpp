#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "blzlog.h"
#include "blzrand.h"
#include "blzmath.h"
#include "blzsim.h"

using namespace std;

double dummyFunction(double d, void* pObject) {
  return d*5. + 6;
}

int main()
{
  BlzLog::setLevel(DEBUG);
  const int N = 16384;
  double r[32678];
  int ISEED1=58;
  int ISEED2=256871;
  BlzLog::debugScalar("ISEED1", ISEED1);
  BlzLog::debugScalar("ISEED2", ISEED2);
  BlzRand::setIX(ISEED1, ISEED2);
  BlzLog::debugScalar("ISEED1", ISEED1);
  BlzLog::debugScalar("ISEED2", ISEED2);
  BlzRand::rnstnr(r, N/2);
  BlzRand::getIX(&ISEED1, &ISEED2);
  BlzLog::debugScalar("ISEED1", ISEED1);
  BlzLog::debugScalar("ISEED2", ISEED2);
  int i;
  for(i=3; i>=0; i--) {
    stringstream ss;
    string msg = "r[";
    int ii = (N/(1<<i))-1;
    ss << ii;
    msg+=ss.str();
    msg+="]";
    //BlzLog::debugScalar(msg, r[ii]);
  }

  // Test the Fourier routine
  const int NN = 8;
  double step = TWOPI/(double)(NN-1);
  double data[2*NN];
  for(i=0; i<NN; i++) {
    // set complex part to zero
    data[i*2+1] = 0.0;
    // set real part to a sine wave
    data[i*2] = sin(i*step);
  }

  //BlzLog::debugComplexVector("data", data, NN);
  BlzMath::fourier(data, NN, 1);
  //BlzLog::debugComplexVector("data", data, NN);

	// DUMMY COMMENT

  // Now test psdsim
  float spsd[BLZSIM_DIM16384]; // where to put the output light curve
  float psdslp = 1.7; // psd slope
  float tinc= 1.93893981; 
  BlzSim* pSim = new BlzSim();
  pSim->psdsim(8192, -psdslp, -psdslp, 1.0, tinc, spsd);
  delete pSim;
  BlzLog::debugScalarPair("spsd[0]", spsd[0], spsd[1]);
  for(i=3; i>=0; i--) {
    stringstream ss;
    string msg = "spsd[";
    int ii = (N/(1<<i))-1;
    ss << ii;
    msg+=ss.str();
    msg+="]";
    if(i!=0)
      BlzLog::debugScalarPair(msg, spsd[ii], spsd[ii+1]);
    else
      BlzLog::debugScalarPair(msg, spsd[ii], spsd[ii-1]);
  }

  BlzSimCommon common;
  double qgresult;
  common.tdust = 1200.;
  common.gammad = 11.473015869669796;
  common.betad = 0.99619423459093981;
  common.freq = 4.65842962e+11;

  // Test qg5() and sdgran()
  qgresult = BlzMath::qg5(0.74151808435155919, 0.99292683212516131, sdgran, &common);
  BlzLog::debugScalar("qgresult", qgresult);
  BlzLog::debugScalar("expected", 1.90360089e-11);

  common.freq = 5.86461544e+11;
  qgresult = BlzMath::qg5(0.74151808435155919, 0.99292683212516131, sdgran, &common);
  BlzLog::debugScalar("qgresult", qgresult);
  BlzLog::debugScalar("expected", 2.97223149e-11);

  // Test seedph()
  double seedphResult;
  double freq = 4.65842962e+11;
  pSim = new BlzSim();
  pSim->common.tdust = 1200.;
  pSim->common.gammad = 11.473015869669796;
  pSim->common.betad = 0.99619423459093981;
  pSim->common.freq = 4.65842962e+11;
  pSim->common.dsnth1 = 0.74151808435155919;
  pSim->common.dsnth2 = 0.99292683212516131;
  seedphResult = pSim->seedph(freq);

  BlzLog::debugScalar("seedphResult", seedphResult);
  BlzLog::debugScalar("expected", 9.51800444e-12);
  delete pSim;

  // Test ajnu()
  // Need to initialize a few arrays first
  static const double gam[] = {1.123009204864502, 1.6150462627410889, 2.6083743572235107, 4.2126455307006836, 6.803617000579834, 10.988156318664551, 17.746381759643555, 28.661226272583008, 46.289207458496094, 74.759208679199219, 120.73958587646484, 195.00001525878906, 196.69918823242188, 198.41317749023438, 200.14210510253906, 201.88607788085938, 203.645263671875, 205.41976928710938, 207.20974731445312, 209.01531982421875, 210.83662414550781, 212.67379760742188, 214.5269775390625, 216.39631652832031, 218.28193664550781, 220.18399047851562, 222.10261535644531, 224.03794860839844, 225.99015808105469, 227.95938110351562, 229.94575500488281, 231.94944763183594, 233.97059631347656, 236.00935363769531, 238.06587219238281, 240.14031982421875, 242.23283386230469, 244.34358215332031, 246.47273254394531, 248.62042236328125, 250.78683471679688, 252.97213745117188, 255.17646789550781, 257.39999389648438};

  static const double edist[] = {521.562622, 252.175659, 96.6792526, 37.0649376, 14.2099771, 5.44782877, 2.08859134, 0.800725341, 0.306982487, 0.117691122, 0.0451204814, 0.0172983147, 0.0163322855, 0.0154056689, 0.0145170139, 0.0136649264, 0.0128480354, 0.0120650455, 0.0113146892, 0.0105957566, 0.00990707427, 0.00924751349, 0.00861598272, 0.00801142678, 0.0074328375, 0.00687923562, 0.00634967862, 0.00584326033, 0.00535909832, 0.00489634927, 0.00445420248, 0.00403186772, 0.00362859014, 0.0032436403, 0.0028763141, 0.00252593099, 0.00219184044, 0.00187341019, 0.0015700314, 0.00128112198, 0.00100611441, 0.000744464516, 0.000495651504, 0.000259169785};

  pSim = new BlzSim();
  double anu = 1e10;
	pSim->common.bperpp = 8.03727818;
  pSim->common.setGgam(gam);
  pSim->common.setEdist(edist);
  double ajnuResult = pSim->ajnu(anu);
  BlzLog::debugScalar("ajnuResult", ajnuResult);
  BlzLog::debugScalar("expected", 26.3672028);
  delete pSim;
}


