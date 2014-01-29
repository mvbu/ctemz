#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <omp.h>
#include <list>
#include "blzutil.h"
#include "blztimer.h"
#include "blzlog.h"
#include "blzrand.h"
#include "blzmath.h"
#include "blzsim.h"
#include "blzsiminputreader.h"

using namespace std;

int main()
{
  BlzLog::setLevel(DEBUG);
  int nThreads = BlzUtil::getNumProcessors();
  const int NUM_THREADS = nThreads;
  const int D68=68;
  const int D44=44;
  double fsscmd_par[D68];
  double fsscmd[D68];
  int i;
  list<int> threadFreqs[NUM_THREADS];

  for(i=0; i<D68; i++) { fsscmd_par[i]=0.0; }
  for(i=0; i<D68; i++) { fsscmd[i]=0.0; }

  // Inputs to ssc() when md=1 (and inu=23 but probably the same for all inu)
  const static double ggam[] = {10, 13.623344859430983, 18.55955251589846, 25.284318436080468, 34.445698959039291, 46.926563584313499, 63.929675877711858, 87.093602123371483, 118.65061767767602, 161.64182824074786, 220.21023698326184, 300, 322.28111280506243, 346.21705223623121, 371.93072288927004, 399.55415753048419, 429.22919504937886, 461.10820876311669, 495.35488881242094, 532.14508266631663, 571.66769805147487, 614.125672942459, 659.73701759351161, 708.73593396243336, 761.37401827449412, 817.92155290121957, 878.66889418749929, 943.92796335313631, 1014.0338481242167, 1089.3465233182408, 1170.2526992177613, 1257.167807223427, 1350.538132982221, 1450.8431079439385, 1558.5977711124262, 1674.355413632009, 1798.7104197883477, 1932.3013190114989, 2075.8140645523858, 2229.985555667809, 2395.6074213994566, 2573.5300853756094, 2764.6671325071957, 2970};
  const static double edist[] =  {6.432732921052045, 3.4659956998109824, 1.8674996053067319, 1.0062201681355147, 0.54215755863379078, 0.29211779657368864, 0.15739484899943296, 0.084805303826482631, 0.045693614516749193, 0.024619997964717713, 0.013265404941001549, 0.0071474810233911334, 0.0056835186043016064, 0.0045165193980876992, 0.0035866196859259862, 0.0028459743135883838, 0.0022563502602039268, 0.0017872011793335794, 0.0014141268945934331, 0.0011176409851825528, 0.0008821849312469314, 0.00069533956607253716, 0.00054719441104330668, 0.00042984333809141967, 0.00033698130373813591, 0.00026358194164099623, 0.00020563983734396801, 0.00015996454026025341, 0.00012401595441285322, 9.5772819712220401e-05, 7.3627652526180297e-05, 5.6302840431783239e-05, 4.27836473336862e-05, 3.2264734431149151e-05, 2.4107482121010329e-05, 1.7805941705657433e-05, 1.2959680846698768e-05, 9.2521347720871591e-06, 6.4333536820803949e-06, 4.3062595161092256e-06, 2.7157033713626476e-06, 1.5397573156969088e-06, 6.8278824451289197e-07, 6.995250169809215e-08};
  const static double snu[] = {10000000000, 15000000000, 43000000000, 86000000000, 100000000000, 230000000000, 316227766016.83795, 562341325190.34912, 1000000000000, 1778279410038.9229, 3162277660168.3794, 5623413251903.4912, 10000000000000, 17782794100389.227, 31622776601683.793, 56234132519034.906, 100000000000000, 177827941003892.28, 316227766016837.94, 562341325190349.06, 1000000000000000, 1778279410038922.8, 3162277660168379.5, 5623413251903491, 10000000000000000, 17782794100389228, 31622776601683792, 56234132519034912, 1e+17, 1.7782794100389229e+17, 3.1622776601683795e+17, 5.6234132519034906e+17, 1e+18, 1.7782794100389228e+18, 3.1622776601683794e+18, 5.6234132519034911e+18, 1e+19, 1.7782794100389229e+19, 3.1622776601683792e+19, 5.6234132519034905e+19, 1e+20, 1.7782794100389229e+20, 3.1622776601683794e+20, 5.6234132519034906e+20, 1e+21, 1.7782794100389228e+21, 3.1622776601683794e+21, 5.6234132519034909e+21, 1e+22, 1.7782794100389228e+22, 3.1622776601683792e+22, 5.6234132519034907e+22, 9.9999999999999992e+22, 1.7782794100389228e+23, 3.1622776601683793e+23, 5.6234132519034907e+23, 9.9999999999999998e+23, 1.7782794100389228e+24, 3.1622776601683796e+24, 5.6234132519034903e+24, 1.0000000000000001e+25, 1.7782794100389228e+25, 3.1622776601683795e+25, 5.6234132519034905e+25, 1e+26, 1.7782794100389227e+26, 3.1622776601683796e+26, 5.623413251903491e+26};
  const static double ssseed[] = {1.7809742638516951e-06, 1.7981302468267946e-06, 1.1007075130410817e-06, 7.1808837363054193e-07, 6.5001489497011419e-07, 3.5546617969438332e-07, 2.7620526374537854e-07, 1.6496973928559461e-07, 9.1597124352101507e-08, 4.6803644500163518e-08, 2.2603926287065257e-08, 1.0380205633740367e-08, 4.3931075537415327e-09, 1.6356033513626085e-09, 4.9625846555116387e-10, 1.075606573313539e-10, 1.3212058056588331e-11, 6.0939402802690962e-13, 5.0425199627599929e-15, 2.067928687338561e-18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  // This looks like the same as snu[] Can probably get rid of this
  const static double nu[] = {10000000000, 15000000000, 43000000000, 86000000000, 100000000000, 230000000000, 316227766016.83795, 562341325190.34912, 1000000000000, 1778279410038.9229, 3162277660168.3794, 5623413251903.4912, 10000000000000, 17782794100389.227, 31622776601683.793, 56234132519034.906, 100000000000000, 177827941003892.28, 316227766016837.94, 562341325190349.06, 1000000000000000, 1778279410038922.8, 3162277660168379.5, 5623413251903491, 10000000000000000, 17782794100389228, 31622776601683792, 56234132519034912, 1e+17, 1.7782794100389229e+17, 3.1622776601683795e+17, 5.6234132519034906e+17, 1e+18, 1.7782794100389228e+18, 3.1622776601683794e+18, 5.6234132519034911e+18, 1e+19, 1.7782794100389229e+19, 3.1622776601683792e+19, 5.6234132519034905e+19, 1e+20, 1.7782794100389229e+20, 3.1622776601683794e+20, 5.6234132519034906e+20, 1e+21, 1.7782794100389228e+21, 3.1622776601683794e+21, 5.6234132519034909e+21, 1e+22, 1.7782794100389228e+22, 3.1622776601683792e+22, 5.6234132519034907e+22, 9.9999999999999992e+22, 1.7782794100389228e+23, 3.1622776601683793e+23, 5.6234132519034907e+23, 9.9999999999999998e+23, 1.7782794100389228e+24, 3.1622776601683796e+24, 5.6234132519034903e+24, 1.0000000000000001e+25, 1.7782794100389228e+25, 3.1622776601683795e+25, 5.6234132519034905e+25, 1e+26, 1.7782794100389227e+26, 3.1622776601683796e+26, 5.623413251903491e+26};
  double restnu = 3162277660168379.5;
  double dopref = 4.6988679517179026;
  double dopref2 = dopref*dopref;
  double dtfact = 0.029775133850520243;

  BlzLog::warnScalar("NUM_THREADS = ", NUM_THREADS);
  // Create objects to time ssc() calls
  BlzTimer sscTimer(NUM_THREADS);
  BlzTimer parallelTimer(1);
  double fq1 = 0.98*nu[6];
  int tid;
  omp_set_num_threads(NUM_THREADS);
  int threadIntervals[NUM_THREADS][2];
  // storage space for which freqs were done by which threads;

  BlzMath::getSubIntervals(7, 68, NUM_THREADS, threadIntervals);
  int inu, md;
  const int MDMAX = 100000;

  //
  // Do in parallel
  //
  BlzIndexTracker *pt = new BlzIndexTracker(7,68);
  BlzSim* pSim = new BlzSim();
  pSim->common.setGgam(ggam);
  pSim->common.setEdist(edist);
  pSim->common.setSnu(snu);
  pSim->common.setSsseed(ssseed);
  pSim->common.cscat = 0.0;
  pSim->common.nuhi = 20; // ssc() now has it's own internal nuhi, so I don't think this matters anymore
  pSim->common.betd = 0.33333333333333331;
  pSim->common.gamd = 1.0606601717798212;
  double fssc1 = pSim->ssc(fq1/dopref)*dopref2/dtfact;
  parallelTimer.start();    

  for(md=0; md<MDMAX; md++) {
    #pragma omp parallel shared(pt), private(tid, inu, restnu)
    {
      tid = omp_get_thread_num();
      // we just use the intervals calculated above to give each thread a different "starting" point
      int previousInu = threadIntervals[tid][0];
      while((inu = pt->getNextIndex(previousInu)) >= 0) {
        restnu = nu[inu-1];
        double anuf = restnu/dopref;
        //cout << "Thread " << tid << " doing frequency " << inu << endl;
        sscTimer.start(tid);
        fsscmd_par[inu-1] = pSim->ssc(anuf)*dopref2/dtfact;
        sscTimer.end(tid);
        //threadFreqs[tid].push_back(inu);
        previousInu = inu;
      }
    } // end pragma parallel

    pt->reset();
  } // end md loop


  parallelTimer.end();

  delete pt;

  /*
  int t;
  list<int>::iterator it;
  for(t = 0; t < NUM_THREADS; t++) {
    cout << "Thread " << t << ": ";
    for(it=threadFreqs[t].begin(); it != threadFreqs[t].end(); ++it) 
      cout << *it << " ";
    cout << endl;
  }
  */

  cout << "% Parallel time: " << parallelTimer.getTotalTime() << " sec" << endl;
  cout << "% Parallel ssc() time: " << sscTimer.getTotalTime() << " sec" << endl;
  cout << "% Parallel ssc() timer detail: " << sscTimer.toString() << endl;
  
  for(i=7; i<=D68; i++) cout << "("<<i<<","<<fsscmd_par[i-1]<<")";
  cout << endl;

  parallelTimer.reset();
  sscTimer.reset();
  delete pSim;

  //
  // Now do serially
  //
  pSim = new BlzSim();
  pSim->common.setGgam(ggam);
  pSim->common.setEdist(edist);
  pSim->common.setSnu(snu);
  pSim->common.setSsseed(ssseed);
  pSim->common.cscat = 0.0;
  pSim->common.nuhi = 20; // ssc() now has it's own internal nuhi, so I don't think this matters anymore
  pSim->common.betd = 0.33333333333333331;
  pSim->common.gamd = 1.0606601717798212;
  fssc1 = pSim->ssc(fq1/dopref)*dopref2/dtfact;
  parallelTimer.start();

  for(md=0; md<MDMAX; md++) {
    for(inu=7; inu<=68; inu++) { // 129  why inu 7?
      double restnu = nu[inu-1];
      //double restnu = nu[23];
      double anuf = restnu/dopref;
      sscTimer.start();
      fsscmd[inu-1] = pSim->ssc(anuf)*dopref2/dtfact;
      sscTimer.end();
    }
  } // end md loop

  parallelTimer.end();

  cout << "% Serial time: " << parallelTimer.getTotalTime() << " sec" << endl;
  cout << "% Serial ssc() time: " << sscTimer.getTotalTime() << " sec" << endl;
  cout << "% Serial ssc() timer detail: " << sscTimer.toString() << endl;

  for(i=7; i<=D68; i++) cout << "("<<i<<","<<fsscmd[i-1]<<")";
  cout << endl;

  cout << endl;

  cout << "Two arrays are ";
  if(!BlzUtil::arraysEqual(fsscmd_par, fsscmd, D68))
    cout << "not ";
  cout << "equal" << endl;

  delete pSim;
}
