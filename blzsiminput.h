#ifndef _INCL_BLZSIMINPUT_H_
#define _INCL_BLZSIMINPUT_H_


class BlzSimInput {
  // This will contain some of the input parameters for a call to BlzSim:run()
 public:
  // These are the input parameters from temz.f and temzinp.txt:
  // read(2,9111)dumdum, zred, dgpc, alpha, p, bave, psdslp,
  //   , uratio, rsize, gmaxmn, gmrat, gmin, betaup, betat, zeta, thlos,
  //   , opang,tdust,dtdist,dtrad,zdist0,vmd
  double nend, zred, dgpc, alpha, p, bave, psdslp;
  double uratio, rsize, gmaxmn, gmrat, gmin, betaup, betat, zeta, thlos;
  double opang, tdust, ldust, dtdist, dtrad, zdist0, vmd;

  // Later we can add other run-time parameters such as total sim time (or number of time steps), etc.
};


#endif // _INCL_BLZSIMINPUT_H_
