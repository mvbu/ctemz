#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "blzlog.h"
#include "blzrand.h"
#include "blzmath.h"
#include "blzsim.h"
#include "blzsiminputreader.h"

using namespace std;

double dummyFunction(double d, void* pObject) {
  return d*5. + 6;
}

int main()
{
  BlzLog::setLevel(DEBUG);
  const int N = 16384;
  double r[32678];
  int ISEED1=58;
  int ISEED2=256871;
  BlzRand randObj;
  BlzLog::debugScalar("ISEED1", ISEED1);
  BlzLog::debugScalar("ISEED2", ISEED2);
  randObj.setIX(ISEED1, ISEED2);
  BlzLog::debugScalar("ISEED1", ISEED1);
  BlzLog::debugScalar("ISEED2", ISEED2);
  randObj.rnstnr(r, N/2);
  randObj.getIX(&ISEED1, &ISEED2);
  BlzLog::debugScalar("ISEED1", ISEED1);
  BlzLog::debugScalar("ISEED2", ISEED2);
  int i;
  for(i=3; i>=0; i--) {
    stringstream ss;
    string msg = "r[";
    int ii = (N/(1<<i))-1;
    ss << ii;
    msg+=ss.str();
    msg+="]";
    //BlzLog::debugScalar(msg, r[ii]);
  }

  // Test BlzRand::rand()
  cout << "BlzRand::rand()" << endl;
  for(i=0; i<44; i++) {
    double r = randObj.rand(0);
    cout << r << ", ";
  }
  cout << endl;
  cout << "BlzRand::rand() test data" << endl;
  randObj.setRandTestMode(true);
  for(i=0; i<100000; i++) {
    double r = randObj.rand(0);
    if(i>99997)  cout << r << ", ";
  }
  cout << endl;
  for(i=0; i<44; i++) {
    double r = randObj.rand(0);
    cout << r << ", ";
  }
  cout << endl;
  randObj.setRandTestMode(false);

  // Test the Fourier routine
  const int NN = 8;
  double step = TWOPI/(double)(NN-1);
  double data[2*NN];
  for(i=0; i<NN; i++) {
    // set complex part to zero
    data[i*2+1] = 0.0;
    // set real part to a sine wave
    data[i*2] = sin(i*step);
  }

  //BlzLog::debugComplexVector("data", data, NN);
  BlzMath::fourier(data, NN, 1);
  //BlzLog::debugComplexVector("data", data, NN);

	// DUMMY COMMENT

  // Now test psdsim
  double spsd[BLZSIM_DIM16384]; // where to put the output light curve
  float psdslp = 1.7; // psd slope
  float tinc= 1.93893981; 
  BlzSim* pSim = new BlzSim();
  pSim->psdsim(16384, -psdslp, -psdslp, 1.0, tinc, spsd);
  delete pSim;
  BlzLog::debugScalarPair("spsd[0]", spsd[0], spsd[1]);
  for(i=3; i>=0; i--) {
    stringstream ss;
    string msg = "spsd[";
    int ii = (N/(1<<i))-1;
    ss << ii;
    msg+=ss.str();
    msg+="]";
    if(i!=0)
      BlzLog::debugScalarPair(msg, spsd[ii], spsd[ii+1]);
    else
      BlzLog::debugScalarPair(msg, spsd[ii], spsd[ii-1]);
  }

  BlzSimCommon common;
  double qgresult;
  common.tdust = 1200.;
  common.gamd = 11.473015869669796;
  common.betd = 0.99619423459093981;
  common.freq = 4.65842962e+11;

  // Test qg5() and sdgran()
  qgresult = BlzMath::qg5(0.74151808435155919, 0.99292683212516131, sdgran, &common);
  BlzLog::debugScalar("qgresult", qgresult);
  BlzLog::debugScalar("expected", 1.90360089e-11);

  common.freq = 5.86461544e+11;
  qgresult = BlzMath::qg5(0.74151808435155919, 0.99292683212516131, sdgran, &common);
  BlzLog::debugScalar("qgresult", qgresult);
  BlzLog::debugScalar("expected", 2.97223149e-11);

  // Test seedph()
  double freq = 159372195086431.97;
  double seedphResult;
  pSim = new BlzSim();
  pSim->common.gamd = 5.4293398294782591;
  pSim->common.betd = 0.98289169615318539;
  pSim->common.tdust = 1200.;
  pSim->common.csth1 = 0.87749686244748315;
  pSim->common.csth2 = 0.63260046099041345;
  seedphResult = pSim->seedph(freq);
  BlzLog::debugScalar("seedphResult", seedphResult);
  BlzLog::debugScalar("expected", 9.51800444e-12);
  delete pSim;

  // Test ajnu()
  // Need to initialize a few arrays first
  static const double gam[] = {1.123009204864502, 1.6150462627410889, 2.6083743572235107, 4.2126455307006836, 6.803617000579834, 10.988156318664551, 17.746381759643555, 28.661226272583008, 46.289207458496094, 74.759208679199219, 120.73958587646484, 195.00001525878906, 196.69918823242188, 198.41317749023438, 200.14210510253906, 201.88607788085938, 203.645263671875, 205.41976928710938, 207.20974731445312, 209.01531982421875, 210.83662414550781, 212.67379760742188, 214.5269775390625, 216.39631652832031, 218.28193664550781, 220.18399047851562, 222.10261535644531, 224.03794860839844, 225.99015808105469, 227.95938110351562, 229.94575500488281, 231.94944763183594, 233.97059631347656, 236.00935363769531, 238.06587219238281, 240.14031982421875, 242.23283386230469, 244.34358215332031, 246.47273254394531, 248.62042236328125, 250.78683471679688, 252.97213745117188, 255.17646789550781, 257.39999389648438};

  static const double edist[] = {521.562622, 252.175659, 96.6792526, 37.0649376, 14.2099771, 5.44782877, 2.08859134, 0.800725341, 0.306982487, 0.117691122, 0.0451204814, 0.0172983147, 0.0163322855, 0.0154056689, 0.0145170139, 0.0136649264, 0.0128480354, 0.0120650455, 0.0113146892, 0.0105957566, 0.00990707427, 0.00924751349, 0.00861598272, 0.00801142678, 0.0074328375, 0.00687923562, 0.00634967862, 0.00584326033, 0.00535909832, 0.00489634927, 0.00445420248, 0.00403186772, 0.00362859014, 0.0032436403, 0.0028763141, 0.00252593099, 0.00219184044, 0.00187341019, 0.0015700314, 0.00128112198, 0.00100611441, 0.000744464516, 0.000495651504, 0.000259169785};

  pSim = new BlzSim();
  double anu = 1e10;
  pSim->common.bperpp = 8.03727818;
  pSim->common.setGgam(gam);
  pSim->common.setEdist(edist);
  double ajnuResult = pSim->ajnu(anu);
  BlzLog::debugScalar("ajnuResult", ajnuResult);
  BlzLog::debugScalar("expected", 26.3672028);
  delete pSim;

  // Test akapnu()
  static const double edist2[] = {207.588089, 100.368889, 38.4794807, 14.7522821, 5.65573835, 2.16830015, 0.831283987, 0.318698138, 0.122182652, 0.0468424521, 0.0179584846, 0.00688493345, 0.00650044205, 0.00613163831, 0.00577794248, 0.0054388023, 0.00511366921, 0.00480203005, 0.00450337958, 0.0042172363, 0.0039431327, 0.00368061941, 0.00342926243, 0.00318864221, 0.0029583571, 0.00273801689, 0.00252724718, 0.00232568663, 0.00213298434, 0.00194880483, 0.00177282514, 0.00160473085, 0.00144422159, 0.00129100704, 0.00114480685, 0.00100535038, 0.000872378412, 0.000745639321, 0.000624891079, 0.000509901613, 0.000400445395, 0.000296305661, 0.000197275149, 0.000103152641};

  pSim = new BlzSim();
  anu = 1e10; // same as before, at the moment 
  pSim->common.bperpp = 12.7397528;
  pSim->common.setGgam(gam);
  pSim->common.setEdist(edist2);
  double akapnuResult = pSim->akapnu(anu);
  BlzLog::debugScalar("akapnuResult", akapnuResult);
  BlzLog::debugScalar("expected", 1.16971385);
  delete pSim;

  // Test ssc()
  static double edist_ssc[] = { 443.1716, 214.273575, 82.1483231, 31.4940643, 12.074213, 4.62901831, 1.77467537, 0.680376053, 0.260842919, 0.10000211, 0.0383388586, 0.0146983732, 0.013877538, 0.0130901923, 0.0123351021, 0.0116110835, 0.0109169716, 0.0102516655, 0.00961408857, 0.00900321081, 0.00841803849, 0.00785760954, 0.00732099777, 0.00680730632, 0.00631567976, 0.00584528456, 0.00539532024, 0.00496501662, 0.00455362396, 0.00416042656, 0.00378473452, 0.00342587661, 0.00308321184, 0.00275612017, 0.00244400324, 0.00214628293, 0.00186240615, 0.00159183599, 0.00133405521, 0.00108856894, 0.000854895043, 0.000632571289, 0.000421154924, 0.000220216491 };
  static double snu[] = { 1e+10, 1.77827942e+10, 3.16227768e+10, 5.62341315e+10, 9.9999998e+10, 1.77827938e+11, 3.16227781e+11, 5.62341347e+11, 9.99999996e+11, 1.77827938e+12, 3.16227768e+12, 5.62341347e+12, 9.99999983e+12, 1.77827941e+13, 3.16227768e+13, 5.62341326e+13, 1e+14, 1.77827937e+14, 3.16227772e+14, 5.6234131e+14, 9.99999987e+14, 1.7782794e+15, 3.16227758e+15, 5.62341303e+15, 1.00000003e+16, 1.77827945e+16, 3.16227769e+16, 5.62341314e+16, 9.99999984e+16, 1.77827933e+17, 3.16227756e+17, 5.62341339e+17, 9.99999984e+17, 1.77827946e+18, 3.16227763e+18, 5.62341319e+18, 9.99999998e+18, 1.77827941e+19, 3.16227768e+19, 5.62341352e+19, 1.00000002e+20, 1.77827939e+20, 3.16227777e+20, 5.62341334e+20, 1.00000002e+21, 1.77827939e+21, 3.16227763e+21, 5.62341348e+21, 9.99999978e+21, 1.77827941e+22, 3.16227769e+22, 5.62341337e+22, 9.99999978e+22, 1.77827946e+23, 3.16227778e+23, 5.62341319e+23, 1.00000001e+24, 1.77827939e+24, 3.16227749e+24, 5.62341326e+24, 9.99999956e+24, 1.77827939e+25, 3.16227761e+25, 5.62341338e+25, 1.00000003e+26, 1.77827941e+26, 3.16227774e+26, 5.62341319e+26 };
    static const double ssseed[]  = { 1.03139553e-06, 4.12530653e-06, 1.72026703e-05, 5.10973186e-05, 5.7645284e-05, 4.3994758e-05, 2.86946997e-05, 1.75751247e-05, 9.08344646e-06, 4.14674287e-06, 1.46940556e-06, 3.13423982e-07, 2.33193926e-08, 3.34458017e-10, 4.20477265e-13, 4.55808116e-18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  pSim = new BlzSim();
  double anuf =  2.86500065e+10;
  pSim->common.nuhi = 16;
  pSim->common.setGgam(gam);
  pSim->common.setEdist(edist_ssc);
  pSim->common.setSnu(snu);
  pSim->common.setSsseed(ssseed);
  double sscResult = pSim->ssc(anuf);
  BlzLog::debugScalar("sscResult", sscResult);
  BlzLog::debugScalar("expected",  0.024821050307809192);
  delete pSim;

  // Test ecdust()
  static const double ggam_ecdust[] = { 0.87952703237533569, 1.6948215961456299, 3.2668538093566895, 6.297025203704834, 12.137832641601562, 23.396284103393555, 45.097515106201172, 86.927726745605469, 167.55757141113281, 322.97564697265625, 622.55181884765625, 1200, 1219.662841796875, 1239.6478271484375, 1259.9603271484375, 1280.6055908203125, 1301.5892333984375, 1322.9166259765625, 1344.5936279296875, 1366.6256103515625, 1389.018798828125, 1411.77880859375, 1434.9117431640625, 1458.423828125, 1482.321044921875, 1506.60986328125, 1531.2967529296875, 1556.3880615234375, 1581.890625, 1607.8109130859375, 1634.156005859375, 1660.932861328125, 1688.1483154296875, 1715.809814453125, 1743.924560546875, 1772.5, 1801.5435791015625, 1831.0631103515625, 1861.0662841796875, 1891.5611572265625, 1922.5556640625, 1954.05810546875, 1986.07666015625, 2018.619873046875 };
   static const double edist_ecdust[] = { 54.3787384, 14.6458302, 3.94187236, 1.06094074, 0.285548359, 0.0768543035, 0.0206850581, 0.00556730898, 0.00149842107, 0.000403294631, 0.00010854529, 2.92145687e-05, 2.6896425e-05, 2.47395292e-05, 2.27332275e-05, 2.0867581e-05, 1.91332601e-05, 1.75215719e-05, 1.60243453e-05, 1.46339853e-05, 1.33433377e-05, 1.21457533e-05, 1.10349883e-05, 1.00052102e-05, 9.05097568e-06, 8.16718057e-06, 7.34905689e-06, 6.59215675e-06, 5.89230831e-06, 5.24562711e-06, 4.64846926e-06, 4.09743552e-06, 3.58935517e-06, 3.12125439e-06, 2.69036309e-06, 2.29409102e-06, 1.93002074e-06, 1.5958916e-06, 1.28959528e-06, 1.00916122e-06, 7.52751987e-07, 5.18650722e-07, 3.05256975e-07, 1.11075792e-07 };
   static const double dustnu[] = { 1.4985991e+11, 1.88662448e+11, 2.37511967e+11, 2.99009901e+11, 3.76431051e+11, 4.73898648e+11, 5.96603109e+11, 7.51078736e+11, 9.45552163e+11, 1.19037965e+12, 1.49859913e+12, 1.88662468e+12, 2.3751198e+12, 2.99009861e+12, 3.76431051e+12, 4.7389868e+12, 5.96603096e+12, 7.51078788e+12, 9.45552255e+12, 1.19037955e+13, 1.49859913e+13, 1.8866252e+13 };
   static const double dusti[] = { 9.51800444e-12, 1.48611574e-11, 2.31122291e-11, 3.57641555e-11, 5.49886907e-11, 8.38589476e-11, 1.26555683e-10, 1.88442997e-10, 2.75783063e-10, 3.94694805e-10, 5.4881516e-10, 7.35197514e-10, 9.38718325e-10, 1.12717435e-09, 1.25207189e-09, 1.26118194e-09, 1.12446863e-09, 8.61036353e-10, 5.45080592e-10, 2.72187495e-10, 1.00743747e-10, 2.54984783e-11 };

  pSim = new BlzSim();
  anuf = 1.70602292e+14;
  pSim->common.setGgam(ggam_ecdust);
  pSim->common.setEdist(edist_ecdust);
  pSim->common.setDustnu(dustnu);
  pSim->common.setDusti(dusti);
  double ecdustResult = pSim->ecdust(anuf);
  BlzLog::debugScalar("ecdustResult", ecdustResult);
  BlzLog::debugScalar("expected", 2.8416587e-7);
  delete pSim;

  // Test polcalc()
  pSim = new BlzSim();
  double b = 0.316276222;
  double bx = 0.150767073;
  double by = 0.240987927;
  double bz = -0.138653368;
  double clos = 0.99098321431122793;
  double slos = 0.13398607746100638;
  pSim->common.bdx = 0.042337138162786593;
  pSim->common.bdy = 0.073330036710039445;
  pSim->common.bdz = 0.97277779711713219;
  pSim->common.gamd = 4.6357107942851155;
  const double expectedChi = 0.94032985;
  double chi = pSim->polcalc(b, bx, by, bz, clos, slos);
  BlzLog::debugScalar("polcalcResult", chi);
  BlzLog::debugScalar("expected", expectedChi);
  delete pSim;

  //
  // Test vdcalc()
  // 
  BlzLog::debug("vdcalc()");
  double vx,vy,vz,sx,sy,sz,vdx,vdy,vdz,vd,gd,eta;
  vx=-0.019240988424703537;
  vy=-0.031763915994522718;
  vz=0.99428764521559587;
  sx=0.086824050218871157;
  sy=0.15038358911704655;
  sz=0.98480777757993587;
  pSim = new BlzSim();
  pSim->vdcalc(vx,vy,vz,sx,sy,sz,&vdx,&vdy,&vdz,&vd,&gd,&eta);
  BlzLog::debugScalarPair("vdx", 0.041984204418010367, vdx);
  BlzLog::debugScalarPair("vdy", 0.073358681218656629, vdy);
  BlzLog::debugScalarPair("vdz", 0.972837232139588060, vdz);
  BlzLog::debugScalarPair("vd", 0.97650215041635646, vd);
  BlzLog::debugScalarPair("gd", 4.6402063567000154, gd);
  BlzLog::debugScalarPair("eta", 2.4875972258196399, eta);
  delete pSim;

  //
  // Test bdcalc()
  //
  /*
  double bxbd, bybd, bzbd,bdxbd,bdybd,bdzbd;
  vx=-0.019240988424703537;
  vy=-0.031763915994522718;
  vz=0.99428764521559587;
  sx=-0.49240406488732119;
  sy=-0.85286842052300182;
  sz=0.17364803833636447;
  bxbd=0.107818834;
  bybd=-0.0286603495;
  bzbd=0.0073269424;
  eta=2.4875972258196399;
  BlzLog::debug("bdcalc()");
  pSim = new BlzSim();
  pSim->bdcalc(vx,vy,vz,sx,sy,sz,bxbd,bybd,bzbd,eta,&bdxbd,&bdybd,&bdzbd);
  BlzLog::debugScalarPair("bdx", 14.07147766, bdxbd);
  BlzLog::debugScalarPair("bdy", -7.848277665, bdybd);
  BlzLog::debugScalarPair("bdz", 0.1261698678, bdzbd);
  delete pSim;
  */

  //
  // Test bcalc()
  //
  /*
  double bparx,bpary,bparz,bprpx,bprpy,bprpz,bpar,bprp;
  vx=-0.019240988424703537;
  vy=-0.031763915994522718;
  vz=0.99428764521559587;
  sx=0.086824050218871157;
  sy=0.15038358911704655;
  sz=0.98480777757993587;
  bx=0.107818834;
  by=-0.0286603495;
  bz=0.0073269424;
  BlzLog::debug("bcalc()");
  pSim = new BlzSim();
  pSim->bcalc(vx,vy,vz,sx,sy,sz,bx,by,bz,
              &bparx,&bpary,&bparz,&bprpx,&bprpy,&bprpz,
              &bpar,&bprp);
  BlzLog::debugScalarPair("bparx", 0.0132327536, bparx);
  BlzLog::debugScalarPair("bpary", 0.0227437261, bpary);
  BlzLog::debugScalarPair("bparz", 0.0134581225, bparz);
  BlzLog::debugScalarPair("bprpx", 0.0945860818, bprpx);
  BlzLog::debugScalarPair("bprpy", -0.0514040738, bprpy);
  BlzLog::debugScalarPair("bprpz", -0.0061311801, bprpz);
  BlzLog::debugScalarPair("bpar", 0.0295550991, bpar);
  BlzLog::debugScalarPair("bprp", 0.107826233, bprp);
  */

  /*  
  //
  // Test run()
  //
  pSim = new BlzSim();
  // Get the input parameters from the input file
  const string inputFile("temzinp.txt");
  BlzSimInput inp;
  BlzSimInputReader inputReader(inputFile);
  inputReader.read(inp);
  // true means run in test mode, second param is number of days to simulate
  pSim->run(inp, 10, true);
  delete pSim;
  */
}


