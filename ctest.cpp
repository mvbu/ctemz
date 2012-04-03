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
  
}


