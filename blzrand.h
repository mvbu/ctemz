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
  static int IX1;
  static int IX2;

 public:
  static void setIX(int ix1, int ix2) { IX1=ix1; IX2=ix2; }
  static void getIX(int* pix1, int* pix2) { *pix1=IX1; *pix2=IX2; }
  static void rnecuy(double u[], int n);
  static void rnstnr(double r[], int n);

 private:
  BlzRand();
  ~BlzRand();
};
    
#endif // _INCL_BLZRAND_H_

