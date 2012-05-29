#ifndef _INCL_BLZRAND_H_
#define _INCL_BLZRAND_H_

#include <string>
#include <sstream>
#include <cstring>

using namespace std;

class BlzRand {
 protected:
  // These are used to store the input/output seed values for each run of 
  // rnecuy(), and therefore of rnstnr() as well
  int IX1;
  int IX2;

 public:
  void setIX(int ix1, int ix2) { IX1=ix1; IX2=ix2; }
  void getIX(int* pix1, int* pix2) { *pix1=IX1; *pix2=IX2; }
  void rnecuy(double u[], int n);
  void rnstnr(double r[], int n);

 public:
  BlzRand();
  ~BlzRand();
};
    
#endif // _INCL_BLZRAND_H_

