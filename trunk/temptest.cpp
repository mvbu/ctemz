

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

main()
{
  int r = rand();
  printf("Hello World %d\n", r);
  printf("RAND_MAX %d\n", RAND_MAX);
  printf("%lf\n", (double)r/RAND_MAX);
	cout << "C++ Hello World" << endl;
}
