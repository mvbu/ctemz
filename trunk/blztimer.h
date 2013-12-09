#ifndef _INCL_BLZTIMER_H_
#define _INCL_BLZTIMER_H_

class BlzTimer {
 public:
  BlzTimer(int nThreads_=1);
  ~BlzTimer();
  static const int MAX_THREADS = 64;
  
  void reset();
  void start(int iThread=0);
  void end(int iThread=0);

  double getThreadTotalTime(int iThread);
  double getTotalTime(); // total of all threads
  
  const char* toString();

 private:
  void init();
  int nThreads;
  double startTimes[MAX_THREADS];
  double endTimes[MAX_THREADS];
  double executionTimes[MAX_THREADS];
  int nCalls[MAX_THREADS];
  double avgTimePerCall[MAX_THREADS];
  
};
    
#endif // _INCL_BLZTIMER_H_
