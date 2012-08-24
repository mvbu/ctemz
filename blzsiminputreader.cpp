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
  filestream >> simInput.nend;
  filestream >> simInput.zred;
  filestream >> simInput.dgpc;
  filestream >> simInput.alpha;
  filestream >> simInput.p;
  filestream >> simInput.bave;
  filestream >> simInput.psdslp;
  filestream >> simInput.uratio;
  filestream >> simInput.rsize;
  filestream >> simInput.gmaxmn;
  filestream >> simInput.gmrat;
  filestream >> simInput.gmin;
  filestream >> simInput.betaup;
  filestream >> dummy; // these dummy reads are to skip past the "d0"
  filestream >> simInput.betat;
  filestream >> dummy;
  filestream >> simInput.zeta;
  filestream >> dummy;
  filestream >> simInput.thlos;
  filestream >> dummy;
  filestream >> simInput.opang;
  filestream >> dummy;
  filestream >> simInput.tdust;
  filestream >> simInput.ldust;
  filestream >> simInput.dtdist;
  filestream >> simInput.dtrad;
  filestream >> simInput.zdist0;
  filestream >> simInput.vmd;
  filestream.close();
}

