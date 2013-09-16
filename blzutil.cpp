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
