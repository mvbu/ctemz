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

int main()
{
  BlzLog::setLevel(DEBUG);
  const int N = 16384;
  double r[32678];
  int ISEED1=58;
  int ISEED2=256871;
  //BlzLog::debugScalar("ISEED1", ISEED1);
  //BlzLog::debugScalar("ISEED2", ISEED2);
  BlzRand::setIX(ISEED1, ISEED2);
  //BlzLog::debugScalar("ISEED1", ISEED1);
  //BlzLog::debugScalar("ISEED2", ISEED2);
  BlzRand::rnstnr(r, N/2);
  BlzRand::getIX(&ISEED1, &ISEED2);
  //BlzLog::debugScalar("ISEED1", ISEED1);
  //BlzLog::debugScalar("ISEED2", ISEED2);
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

  // Now test psdsim
  float spsd[BLZSIM_DIM16384]; // where to put the output light curve
  float psdslp = 1.7; // psd slope
  float tinc= 1.93893981; 
  BlzSim::psdsim(BLZSIM_DIM16384, -psdslp, -psdslp, 1.0, tinc, spsd);
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
}


