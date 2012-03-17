
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "blzlog.h"
#include "blzrand.h"

using namespace std;

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
		cout << ii << " ";
		ss << ii;
		msg+=ss.str();
		msg+="]";
		BlzLog::debugScalar(msg, r[ii]);
	}
}


