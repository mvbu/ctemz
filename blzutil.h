#ifndef _INCL_BLZUTIL_H_
#define _INCL_BLZUTIL_H_

class BlzUtil {
 public:
  static int getNumProcessors();
  template<class T> 
    static bool arraysEqual(T array1[], T array2[], int nSize) {
    int i;
    for(i=0; i< nSize; i++) {
      if(array1[i] != array2[i])
        return false;
    }
    return true;
  }

 private:
  BlzUtil();
  ~BlzUtil();
};

// Class to assign indices to the caller and keep track of which indices are taken
class BlzIndexTracker {
 public:
  BlzIndexTracker(int startIndex, int endIndex);
  ~BlzIndexTracker();
  void reset();

  // Get next unused index. Pass in current index in case it's useful.
  // Probably want to try and give the next consecutive index if it's available
  int getNextIndex(int currentIndex=-1);

 private:
  void init();
  bool *pIndices;
  int range, iStart, iEnd;
};

#endif // _INCL_BLZUTIL_H_

