#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include "blzutil.h"

using namespace std;

int BlzUtil::getNumProcessors()
{
  // This will not work on OS X (the file isn't there)
  int retVal = 1;
  try {
    FILE * fp;
    char res[128];
    fp = popen("/bin/cat /proc/cpuinfo |grep -c '^processor'","r");
    fread(res, 1, sizeof(res)-1, fp);
    fclose(fp);
    retVal = atoi(res);
  }
  catch(int n) {
    // just return the default value
  }

  return retVal;
}

BlzIndexTracker::BlzIndexTracker(int startIndex, int endIndex):
  iStart(startIndex),iEnd(endIndex),range(endIndex-startIndex+1)
{
  // allocate enough space for the indices
  pIndices = (bool*)malloc(range*sizeof(bool));
  init();
}

BlzIndexTracker::~BlzIndexTracker() 
{
  if(pIndices != NULL) 
    free(pIndices);
}


void BlzIndexTracker::init()
{
  int i;
  for(i=0; i<range; i++) {
    pIndices[i] = false;
  }
}

void BlzIndexTracker::reset() { init(); }

int BlzIndexTracker::getNextIndex(int currentIndex)
{
  if((currentIndex >= iStart) && (currentIndex < iEnd)) {
    // If a valid index was passed, try to return the next consecutive index.
    // if it is "taken", then just return the lowest index that is not already taken.
    if(!pIndices[currentIndex-iStart+1]) {
      pIndices[currentIndex-iStart+1] = true;
      return currentIndex+1;
    }
  }

  // Just try to return the lowest index that is not already taken
  int i;
  for(i=iStart; i<=iEnd; i++) {
    if(!pIndices[i-iStart]) {
      pIndices[i-iStart] = true;
      return i;
    }
  }

  // If all indices are taken, return -1
  return -1;
}
