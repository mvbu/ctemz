#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "blzlog.h"
#include "blzmath.h"
#include "blzutil.h"

using namespace std;

int main()
{
  
  double array1[] = { 1.0, 3.0, 6.2, 1.4, 3.3};
  double array2[] = { 1.0, 3.0, 6.2, 1.4, 3.3};
  double array3[] = { 1.0, 3.0, 6.2, 1.5, 3.3};

  printf("Array 1 equal to Array 2: %d\n", BlzUtil::arraysEqual(array1, array2, 5));
  printf("Array 2 equal to Array 3: %d\n", BlzUtil::arraysEqual(array2, array3, 5));
  printf("Array 1 equal to Array 3: %d\n", BlzUtil::arraysEqual(array1, array3, 5));
  printf("Expect: true, false, false\n");

  const int START = 4;
  const int END = 9;

  BlzIndexTracker *pt = new BlzIndexTracker(START, END);

  int i, count = END - START;

  for(i=0; i<count+1; i++) {
    int n = pt->getNextIndex(5);
    printf("i= %d, n = %d\n", i, n);
  }

  delete pt;

}
