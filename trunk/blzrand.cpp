#include <iostream>
#include <cmath>
#include "blzlog.h"
#include "blzrand.h"

int BlzRand::IX1 = 123456;
int BlzRand::IX2 = 654321;

void BlzRand::rnstnr(double r[], int n)
{
  int i=0;
  const int UDIM = 2;
  double u[UDIM];

  do {
    rnecuy(u, UDIM);
    double V1=2.*u[0]-1.;
    double V2=2.*u[1]-1.;
    double S=V1*V1+V2*V2;
    if(S >= 1.)
      continue;
    double ROOT=sqrt(-2.*log(S)/S);
    r[i]=V1*ROOT;
    if(i < n-1)
      r[i+1]=V2*ROOT;
    i=i+2;
  } while(i <= n-1);
}

void BlzRand::rnecuy(double u[], int n) 
{
  const int M1=2147483563,M1MIN1=2147483562,IA1=40014, IQ1=53668,IR1=12211;
  const int M2=2147483399,IA2=40692,IQ2=52774,IR2=2791;
  const double XM1INV=4.656613E-10;
  BlzLog::psychoScalar("before IX1", IX1);
  BlzLog::psychoScalar("before IX2", IX2);
  int i;
  for(i=0; i<n; i++) {
    // Produce integer random number X1 from first MLCG
    int k = IX1/IQ1;
    IX1=IA1*(IX1-k*IQ1)-k*IR1;
    if (IX1 < 0) 
      IX1=IX1+M1;
    // Produce integer random number X2 from second MLCG
    k=IX2/IQ2;
    IX2=IA2*(IX2-k*IQ2)-k*IR2;
    if (IX2 < 0) 
      IX2=IX2+M2;
    // Combine
    int IZ=IX1-IX2;
    if(IZ < 1) 
      IZ=IZ+M1MIN1;
    // Normalize and transform to floating point
    u[i]=double(IZ)*XM1INV;
  }

  BlzLog::psychoScalar("after IX1", IX1);
  BlzLog::psychoScalar("after IX2", IX2);
}
