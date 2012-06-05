#ifndef _INCL_BLZSIMINPUT_H_
#define _INCL_BLZSIMINPUT_H_


class BlzSimInput {
  // This will contain some of the input parameters for a call to BlzSim:run()
 public:
  double zred, dgpc, alpha, p,bavg, neg_psdslope, e_e_over_e_b, rsie, gmaxmn, gmrat, gmin, betaup, betat;
  double zeta, thlos, opang,dustTemp,dustTorusDist, torusXsectionRadius, zdist0, vmd;
};


#endif // _INCL_BLZSIMINPUT_H_
