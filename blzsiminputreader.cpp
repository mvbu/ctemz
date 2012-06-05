#include <iostream>
#include <fstream>
#include "blzsiminputreader.h"

using namespace std;

BlzSimInputReader::BlzSimInputReader(const std::string _filepath)
  : filepath(std::string(_filepath))
{
}

void BlzSimInputReader::read(BlzSimInput& simInput)
{
  ifstream filestream(filepath.c_str());
  if(!filestream.is_open()) {
    throw new exception();
  }

  // dummy buffer
  char dummy[32768];
  // read in and discard the first line (just some descriptive text)
  filestream.getline(dummy, 32768);
  cout << "Read in: " << dummy << endl;
  filestream >> simInput.zred;
  filestream >> simInput.dgpc;
  filestream >> simInput.alpha;
  filestream >> simInput.p;
  filestream >> simInput.bavg;
  filestream >> simInput.neg_psdslope;
  filestream >> simInput.e_e_over_e_b;
  filestream >> simInput.rsie;
  filestream >> simInput.gmaxmn;
  filestream >> simInput.gmrat;
  filestream >> simInput.gmin;
  filestream >> simInput.betaup;
  filestream >> dummy;
  filestream >> simInput.betat;
  filestream >> dummy;
  filestream >> simInput.zeta;
  filestream >> dummy;
  filestream >> simInput.thlos;
  filestream >> dummy;
  filestream >> simInput.opang;
  filestream >> dummy;
  filestream >> simInput.dustTemp;
  filestream >> simInput.dustTorusDist;
  filestream >> simInput.torusXsectionRadius;
  filestream >> simInput.zdist0;
  filestream >> simInput.vmd;
  filestream.close();
}

