#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sys/time.h>
#include "blztimer.h"

using namespace std;

BlzTimer::BlzTimer(int nThreads_): nThreads(nThreads_) {
  if(nThreads>MAX_THREADS)
    throw new std::exception();
  init();
}

BlzTimer::~BlzTimer() {}

void BlzTimer::init() {
  int i;
  for(i=0;i<MAX_THREADS;i++) {
    startTimes[i] = 0.0;
    endTimes[i] = 0.0;
    executionTimes[i] = 0.0;
    nCalls[i] = 0;
    avgTimePerCall[i] = 0.0;
  }
}

void BlzTimer::reset() { init(); }

void BlzTimer::start(int iThread) {
  if(iThread >= nThreads)
    throw new std::exception();

  timeval tv;
  gettimeofday(&tv, NULL);
  startTimes[iThread] = tv.tv_sec + (1e-6 * tv.tv_usec);
  //cout << startTimes[iThread] << endl;
}

void BlzTimer::end(int iThread) {
  if(iThread >= nThreads)
    throw new std::exception();

  timeval tv;
  gettimeofday(&tv, NULL);
  endTimes[iThread] = tv.tv_sec + (1e-6 * tv.tv_usec);
  double executionTime = (endTimes[iThread] - startTimes[iThread]);
  //cout << "executionTime=" << executionTime << endl;
  executionTimes[iThread] += executionTime; // add onto whatever was already there
  nCalls[iThread]++;
  avgTimePerCall[iThread] = executionTimes[iThread]/static_cast<double>(nCalls[iThread]);
  //cout << startTimes[iThread] << " " << endTimes[iThread] << " " << executionTimes[iThread] << endl;
}

double BlzTimer::getThreadTotalTime(int iThread) {
  return executionTimes[iThread];
}

double BlzTimer::getTotalTime() {
  int i;
  double totalTime = 0.0;
  for(i=0;i<nThreads;i++) {
    totalTime += getThreadTotalTime(i);
  }
  return totalTime;
}

const char* BlzTimer::toString()
{
  string s("BlzTimer ");
  char buf[64];

  int i;
  for(i=0; i<nThreads; i++) {
    sprintf(buf, "[%d,", i);
    s.append(buf);
    sprintf(buf, "%lfms,", executionTimes[i]*1000.);
    s.append(buf);
    sprintf(buf, "(%d,%lfms)]", nCalls[i], avgTimePerCall[i]*1000);
    s.append(buf);
  }
  sprintf(buf, "[Total: %lfms]", getTotalTime()*1000.);
  s.append(buf);
  return s.c_str();
}
