c     Turbulent Extreme Multi-Zone (TEMZ) model for blazar multi-waveband
c       light curves, SEDs, and (synchrotron) polarization 
c     Program that calculates the light curve and dynamic SED of a jet with turbulent
c       plasma passing through a conical standing shock
c     This version includes variations on times as short as the time for a cell to
c       cross the shock, as well as longer-term variations according to a given PSD slope
c     Includes synchrotron radiation and IC of seed photons from a dust torus
c       and from a Mach disk
      dimension fixedRandData(100000)
      dimension itarra(3),
     , pq(400,1140,68),pu(400,1140,68),fpol(400,1140,68),
     , flux(400,1141,68),gammin(400,1141),gammax(400,1141),
     , xcell(1141),gmax0(400,1141),egam(400,1141,44),
     , phcell(1141),rcell(1141),ycell(1141),cosph(1141),
     , sinph(1141),bperp(1141),bfield(1141),n0(400,1141),
     , zcell(400,1141),betadx(400,1141),betady(400,1141),
     , betadz(400,1141),betad(400,1141),gammad(400,1141),
     , betaux(1141),betauy(1141),betauz(1141),
     , betau(1141),gammau(1141),spsdx(131072),spsd(131072),
     , nu(68),bx(400,1141),by(400,1141),bz(400,1141),
     , delta(400,1141),enofe(400,1141,44),edist(44),fsync(68),
     , ggam(44),dustnu(22),dusti(22),ididg(1141),
     , gcnt(44),igcnt(44),fsynmd(68,130000),nouter(1140),
     , fsscmd(68,131000),fmdall(68),deltmd(1140),dmd(1140),
     , tlf(200),betamx(1140),betamy(1140),
     , betamz(1140),betamr(1140),gamamr(1140),bmdx(130000),
     , bmdy(130000),bmdz(130000),bmdtot(130000),tlf1(1140),
     , flsync(400,1141,68),flcomp(400,1141,68),absorb(68,130000),
     , fsync2(68),alphmd(68,130000),dustii(451,22),dustf(451,22),
     , alfmdc(68,130000),syseed(68),scseed(68),snu(68),ssseed(68),
     , flec(400,1141,68),flssc(400,1141,68),mdd(130000),useed(451),
     , phots(68),phalph(68),icelmx(1141),imax(1141),seedpk(451),
     , abexmd(68,130000),psi(1140),sinpsi(1140),cospsi(1140),
     , tanpsi(1140),tauxmd(68),phi1(1141),phi2(1141),
     , theta1(1141),theta2(1141),ididb(1141),bfrac(1141),
     , ididbd(1141),bfracd(1141),phi2d(1141),thet2d(1141),
     , phi1d(1141),thet1d(1141),angrot(1141),cosrot(1141),
     , bu1x(1141),bu1y(1141),bu1z(1141),bu2x(1141),bu2y(1141),
     , bu2z(1141),cu1x(1141),cu1y(1141),cu1z(1141),angrtd(1141),
     , cosrtd(1141),bu1xd(1141),bu1yd(1141),bu1zd(1141),bu2xd(1141),
     , bu2yd(1141),bu2zd(1141),cu1xd(1141),cu1yd(1141),cu1zd(1141),
     , tdelr(451),zsmd(130000)
      character*10 dumdum, fstat
      character*15 filnam
      real*8 pol,pqcum,pucum,pq,pu,pmean,polc,ai2,
     ,  pcum,tanv0,cosv0,betaup,gamup,beta,sinz,cosz,phcell,
     ,  cosph,sinph,thlos,betau,gammau,opang,tanop,cosop,sinop,
     ,  zeta,tanz,betad,betadx,betady,betadz,slos,clos,
     ,  eta,tanxi,xi,betacs,gammad,bdx,bdy,bdz,
     ,  dcsth1,dcsth2,dth1,dth2,dsnth1,dsnth2,n0ave,
     ,  betamd,betarl,betamx,betamy,betamz,betamr,gamamr,
     ,  cosmd,tanmd,bup,gmax0,ustob,tlfact,delt,
     ,  gamb,glow,gmrat,gmratl,gmratm,ggam,tloss,t1,t2,tlmin,
     ,  eterm1,eterm2,gammax,glim,t2max,phots,betaux,betauy,
     ,  betauz,psi,cospsi,sinpsi,betat,
     ,  tanpsi,betd,gamd,bd2,gm1,angm,bmparx,bmpary,bmparz,
     ,  bmprpx,bmprpy,bmprpz,dtnth1,dtnth2,dotprd,
     ,  betupx,betupy,betupz,betatx,betaty,betatz,thetat,
     ,  sintht,costht,phit,btprpx,btprpy,btprpz,
     ,  btparx,btpary,btparz,sx,sy,sz,slx,sly,slz,
     ,  anx,any,anz,betadd,gammdd,csth1,csth2
      real*4 n0mean,n0,nu,ldust
      integer dstart
      real fixedRandData
      integer fixedRandFlag, fixedRandFileOpened, fixedRandCounter
      real*4 daysToSimulate
      character*5 :: dpath='maps/'
c      common/ci/i,j
      common/cvel/bdx,bdy,bdz,gamd,betd
      common/cparm/zred1,bfld,bperpp
      common/cdust/csth1,csth2,dcsth1,dcsth2,dsang,tdust
      common/cssc/snu,ssseed,nuhi,cscat
      common/cdist/ggam,edist
      common/cseed/dustnu,dusti,csang
      common/cang/cosz,sinz
      common/crite/iwrite
      common/cinput/daysToSimulate
      ! Added by me (MSV, June 2012) for test mode ("generates" the same sequence of rand numbers every time)
      common/cfixedrand/fixedRandFlag, fixedRandFileOpened, fixedRandData, fixedRandCounter
      nTestOut = 0
      ! Parse the command line options and their arguments
      daysToSimulate = 0.4 ! Default value
      fixedRandFlag = 0

      ! Replace default value(s) above with values specified on command line
      dummy = parseArgs(0)

      fstat = 'new'

      if(fixedRandFlag.ne.0) then
         fstat = 'replace'
      end if

      open (2,iostat=ios, err=9000, file='temzinp.txt',
     ,   status='old')
      open (3,iostat=ios, err=9000, file='temzspec.txt',
     , status=fstat)
      open (4,iostat=ios, err=9000, file='temzlc.txt',
     , status=fstat)
      open (5,iostat=ios, err=9000, file='temzpol.txt',
     , status=fstat)
      open (6,iostat=ios, err=9000, file='temzcheck.txt',
     , status=fstat)
      it=0
      if(nTestOut.ne.0) open (9,iostat=ios, err=9000, file='testout.txt',
     , status=fstat)
c     Input file format: 1st line characters, then values of parameters,
c       one per line
      read(2,9111)dumdum, nend, zred, dgpc, alpha, p, bave, psdslp,
     , uratio, rsize, gmaxmn, gmrat, gmin, betaup, betat, zeta, thlos,
     , opang,tdust,ldust,dtdist,dtrad,zdist0,vmd
      close(2)
      call itime(itarra)
      ita=itarra(1)+itarra(2)+itarra(3)
      iseed = randproto(ita)
      randstart=randproto(iseed)
c     icells along axial direction, jcells along transverse direction
c      icells=200
      icells=400 ! Can lower this for testing, but Marscher currently has this at 400
      jcells=3*nend*(nend-1)+1
      ancol=2*nend-1
      rbound=ancol*rsize
c     zrat later defines zsize = zrat*rsize/tanz
c       If zrat is decrease, size of dust arrays needs to be increase above 451
      zrat=0.2
      nzid=nend/zrat
      amppsd=6.0
      ndim=131072
c      ndim=65536
      andim=ndim
c      mdmd=17000
c      mdmax=55000
      mdmd=41820
      mdmax=120000
c     An SED will be printed out every ispecs time steps
      ispecs=1
      ispec=1
c     Set up frequencies
      do 1 inu=1,68
      phots(inu)=0.0
      nu(inu)=10.0**(10.0+0.25*(inu-1))
      if(inu.eq.2)nu(inu)=1.5e10
      if(inu.eq.3)nu(inu)=4.3e10
      if(inu.eq.4)nu(inu)=8.6e10
    1 if(inu.eq.6)nu(inu)=2.3e11
c     alpha = (s-1)/2 is the underlying spectral index, where s = slope of 
c       electron E dist.
c     gmaxmn,gmaxmx are the min/max values of initial gamma_max of
c       electrons in a cell
      gmaxmx=gmrat*gmaxmn
c     gmin is the minimum value of initial gamma of electrons
c     2p is the slope of the volume vs. initial gamma-max law over all cells,
c          V = V0*gamma_max**(-2p)
      pexp=-2.0*p
c     Set up compilation of distribution of initial gamma_max values
      gmrf=gmaxmx/gmaxmn
      gmrfl=alog10(gmrf)/43.0
      do 2 ie=1,44
      gfl=alog10(gmaxmn)+gmrfl*(ie-1)
      gcnt(ie)=10.0**gfl
      igcnt(ie)=0
    2 continue
      yr=3.16e7
      rad=180.0/3.1415926
      pi=3.14159265358979
      pio2=1.57079632679490
      sq2=1.41421
      sq3=1.73205
      parsec=3.086e18
      c=3.0e10
      emc2=8.186e-7
      cc2=1.29e-9
      SMALL_FMDALL=1e-30
c     fgeom is multiplier to take into account that grid of cells occupies less
c       than a full cross-section of radius (2*nend-1)*rcell
      fgeom=1.33
      emfact=3.74e-23
      c6gam=8.2e-21
      amjy=1.0e-26
      zred1=1.0+zred
      sen=2.0*alpha+1.0
      zeta=zeta/rad
      opang=opang/rad
c     uratio is the user-specified ratio of energy density of electrons 
c        to that of mag. field
c     Set the normalization of the electron energy distribution accordingly:
      n0ave=uratio*bave**2*(sen-2.0)/(8.0*pi*emc2)/
     , (gmin**(2.0-sen)-gmaxmn**(2.0-sen))
c     Line of sight
      thlos=thlos/rad
      slos=dsin(thlos)
      clos=dcos(thlos)
      sinz=dsin(zeta)
      cosz=dcos(zeta)
      tanz=sinz/cosz
      gamup=1.0d0/dsqrt(1.0d0-betaup**2)
c     Compression ratio of Mach disk
c       Ultra-relativistic eq. of state assumed, so compression ratio is that
c       given by Hughes, Aller, & Aller (1989, ApJ, 341, 54)
      etac=dsqrt(8.0d0*gamup**4-1.70d1*gamup**2+9.0d0)/gamup
c     Speed downstream of conical shock if turbulent velocity is ignored;
c        for setting cell length zsize and for first estimate of time delay
c      betadd=dsqrt((1.0d0-(betaup*cosz)**2)**2+9.0d0*
c     ,(betaup**2*cosz*sinz)**2)/(3.0*betaup*sinz)
      sx=-sinz
      sy=0.0
      sz=cosz
      call vdcalc(0.0d0,0.0d0,betaup,sx,sy,sz,anx,any,anz,
     ,  betadd,gammdd,eta)
      if(betadd.gt.0.57735.and.betadd.lt.betaup)go to 7891
      write(6,9891)betaup,betadd,sinz,cosz
      go to 9000
c 7891 gammdd=1.0d0/dsqrt(1.0d0-betadd**2)
 7891 betd=betadd
      gamd=gammdd
c     Length of a cylindrical cell in pc
      zsize=0.2*rsize/tanz ! 0.2 can be changed to 2.0 for faster runs/testing
      volc=pi*rsize**2*zsize
c     Length and volume of cell in plasma proper frame
      zsizep=zsize/gamd
      volcp=volc/gamd
      svmd=sqrt(vmd)
      delobs=1.0d0/(gamd*(1.0d0-betd*clos))
      dstart=mdmd+icells*delobs
c     Time step in observer's frame in days
      dtfact=(1.0d0-betd*clos)/betd
      dtime=1190.0*zsize*dtfact*zred1
      itlast=daysToSimulate/dtime
      mdrang=0.5*(1.0/dtfact+1.0)
c     Distance of shock from axis and apex of conical jet
      tanop=dtan(opang)
c     Next line is specific to the selected number of cells per slice
      rshock=(2*nend-1)*rsize
c     Distance of Mach disk from z value where conical shock intersects jet boundary
      zshock=rshock/tanz
c     Distance of Mach disk from from vertex of jet cone
      zsvtex=rshock/tanop
c      write(3,8888)iseed,alpha
      expon=alpha+1.0
      exp1=0.5*expon-1.0
      anorm=(1.0+alpha)/(alpha+5.0/3.0)
c     Computation of IR seed photon density from dust torus as function of distance down jet
c     Area filling factor of hot dust, from IR luminosity. 1.05E8 = 1E45/(parsec**2)
      filld=1.05e8*ldust/(5.67e-5*2.0*pi**2*dtrad*dtdist*tdust**4)
      do 333 id=1,(icells+nend/zrat)
      zdist=zdist0+(id-nend/zrat)*zsize+zshock
c     Calculate min & max angles of dust torus in plasma frame
      dphi1=asin(dtrad/sqrt(zdist**2+dtdist**2))
      dphi2=atan(dtdist/zdist)
      dth1=dphi2-dphi1
      dth2=dphi2+dphi1
      csth1=dcos(dth1)
      csth2=dcos(dth2)
      dcsth1=-(csth1-betadd)/(1.0d0-csth1*betadd)
      dcsth2=-(csth2-betadd)/(1.0d0-csth2*betadd)
c     Doppler factor of dust torus emission in frame of cells, used
c       to estimate frequency of peak intensity in plasma frame
      tdel=gammdd*(1.0d0-betadd*0.5*(csth1+csth2))
      tdelr(id)=tdel
      dsang=dth2-dth1
c     Calculate seed photon field from dust emission in plasma frame
c     Peak frequency of dust thermal emission for part of torus closest to shock
      seedpk(id)=5.88e10*tdust*tdel
c     Use this to set the frequency array of the seed photons
      do 3 inu=1,22
      dustnu(inu)=seedpk(id)*10**(-1.4+(inu-1)*0.1)
      dusti(inu)=filld*seedph(dustnu(inu))
      if(dusti(inu).lt.1.0e-20)dusti(inu)=1.0e-20
      dustf(id,inu)=dustnu(inu)
    3 dustii(id,inu)=dusti(inu)
c      do 3331 inu=1,22
c 3331 write(6,9331)id,seedpk(id),dustnu(inu),dusti(inu)
c 9331 format(i6,1p3e12.3)
  333 continue
      do 336 id=1,icells+nzid
      dflux=0.0
      do 6 inu=2,22
      if(dustii(id,inu-1).le.0.0.or.dustii(id,inu).le.0.0)go to 4
      a=alog10(dustii(id,inu)/dustii(id,inu-1))/
     ,  alog10(dustnu(inu)/dustnu(inu-1))
      dflux=dflux+dustii(id,inu-1)/(a+1.0)*dustnu(inu-1)*
     ,  ((dustnu(inu)/dustnu(inu-1))**(a+1.0)-1.0)
      go to 5
    4 dflux=dflux+0.5*(dustii(id,inu)+dustii(id,inu-1))*
     ,  (dustnu(inu)-dustnu(inu-1))
    5 continue
    6 continue
      useed(id)=4.0*pi*dflux/c
c      write(6,9994)id,nzid,useed(id)
  336 continue
      gmratl=log10(gmaxmx/gmin)/40.0
cccccc   TEST
c      do 824 ie=1,44
c      ggam(ie)=gmin*10.0**(gmratl*(ie-4))
c      edist(ie)=ggam(ie)**(-sen)
c  824 continue
c      bperpp=1.0
c      testj=ajnu(1.0e12)
c      write(5,9996)testj
cccccc  END TEST
      write(3,6665)nend,zred,dgpc,alpha,p,bave,uratio,rsize,gmaxmn,
     , gmrat,gmin,betaup,betat,(zeta*rad),(thlos*rad),(opang*rad),
     , tdust,ldust,dtdist,dtrad,zdist0,useed(5),psdslp,vmd,filld
      write(4,6665)nend,zred,dgpc,alpha,p,bave,uratio,rsize,gmaxmn,
     , gmrat,gmin,betaup,betat,(zeta*rad),(thlos*rad),(opang*rad),
     , tdust,ldust,dtdist,dtrad,zdist0,useed(5),psdslp,vmd,filld
      write(5,6665)nend,zred,dgpc,alpha,p,bave,uratio,rsize,gmaxmn,
     , gmrat,gmin,betaup,betat,(zeta*rad),(thlos*rad),(opang*rad),
     , tdust,ldust,dtdist,dtrad,zdist0,useed(5),psdslp,vmd,filld
 6665 format('#No. of cells on each side of hexagonal grid: ',i3/
     ,  '#redshift: ',f5.3,2x,'Distance in Gpc: ',f5.3/
     ,  '#spectral index: ',f4.2,2x,'filling factor exponent: ',f4.2/
     ,  '#mean unshocked magnetic field: ',f5.3,2x,'ratio of electron ',
     ,  'to mag. energy: ',e8.2/'#cell radius (pc): ',f6.3/
     ,  '#Min. value of gamma_max: ',f8.1,2x,'ratio of max. to min. ',
     ,  'values of gamma_max: ',f5.1/'#gamma_min: ',f6.1,/
     ,  '#upstream laminar velocity: ',f9.5,'c',2x,
     ,  'upstream turbulent velocity: ',f9.5,'c'/'#shock angle: ',f6.3,
     ,  2x,'viewing angle: ',f6.3,2x,'opening angle: ',f6.3/
     ,  '#Dust temperature: ',f6.0,2x,'dust luminosity ',f5.2,
     ,  'x10**45 erg/s'/'#distance of center of dust ',
     ,  'torus from black hole: ',f3.1,' pc'/'#radius of torus: ',
     ,  f3.1,' pc   ','Distance of shock from central ',
     ,  'engine: ',f5.2,' pc'/'#Energy density of seed photons in ',
     ,  'plasma frame: ',e9.2/'#-Slope of PSD: ',f5.1,
     ,  '     Area of Mach disk relative to other cell',1pe9.2/
     ,  '# Area filling factor of dust emission: ',e9.2)
c     Set up variation of energy density according to PSD; see Done et al.,
c       1992, ApJ, 400, 138, eq. B1
      tinc=dtime
      call psdsim(ndim,-psdslp,-psdslp,1.0,tinc,spsdx)
      if(nTestOut.eq.6) then
      do ni = 1,ndim
         write(9,14007) ni,spsdx(ni)
      end do
      endif
      psdsum=0.0
      psdsig=0.0
      ipulse=dstart+100+ip0
      spexp=1.0/(0.5*expon+1.0)
      if(nTestOut.eq.6) then
         write(9,14009) 'spexp',spexp
      end if
      do 4997 ip=1,ndim
c     Normalize spsd by standard deviation
      psdsig=psdsig+spsdx(ip)*spsdx(ip)/(andim-1.0)
 4997 continue
      if(nTestOut.eq.6) then
         write(9,14009) 'psdsig^2',psdsig
      end if
      psdsig=sqrt(psdsig)
      if(nTestOut.eq.6) then
         write(9,14009) 'psdsig',psdsig
      end if
      do 4998 ip=1,ndim
cccc  Next section is only for testing purposes
cccc   If using it, comment out call psdsim() above
c      spsdx(ip)=1.0
c  Add pulse of high energy density
c      expsd=3.0-(ipulse-10-ip)/4.0
c      if(ip.gt.(ipulse-10))expsd=3.0+(ipulse-10-ip)/4.0
c      if(ip.ge.(ipulse-20).and.ip.le.(ipulse))spsdx(ip)=
c     ,  exp(0.4*expsd)*spsdx(ip)
c      expsd=3.0-(ipulse+5-ip)/1.5
c      if(ip.gt.(ipulse+5))expsd=3.0+(ipulse+5-ip)/1.5
c      if(ip.ge.(ipulse).and.ip.le.(ipulse+10))spsdx(ip)=
c     ,  exp(1.5*expsd)*spsdx(ip)
c      expsd=3.0-(ipulse+60-ip)/4.0
c      if(ip.gt.(ipulse+60))expsd=3.0+(ipulse+60-ip)/4.0
c      if(ip.ge.(ipulse+50).and.ip.le.(ipulse+70))spsdx(ip)=
c     ,  exp(0.4*expsd)*spsdx(ip)
c      expsd=3.0-(ipulse+95-ip)/4.0
c      if(ip.gt.ipulse+95)expsd=3.0+(ipulse+95-ip)/4.0
c      if(ip.ge.(ipulse+85).and.ip.le.(ipulse+105))spsdx(ip)=
c     ,  exp(0.35*expsd)*spsdx(ip)
cccc
c     Normalize spsd by standard deviation and take exponential of result
c       to get amplitude of flux variation
      spsd(ip)=exp(amppsd*spsdx(ip)/psdsig)
      if(nTestOut.eq.6) then
         write(9,14007) ip,spsd(ip)
      end if
c     Need to scale n0 and B by a different factor to get the desired amplitude
      spsd(ip)=spsd(ip)**(spexp)
      if(nTestOut.eq.6) then
         write(9,14007) ip,spsd(ip)
      end if
c     Average of 10 time steps to smooth variations so that discreteness
c       of columns of cells does not cause artificial spikes of flux
      if(ip.gt.9)spsd(ip)=0.1*(spsd(ip)+spsd(ip-1)+spsd(ip-2)+
     ,  spsd(ip-3)+spsd(ip-4)+spsd(ip-5)+spsd(ip-6)+
     ,  spsd(ip-7)+spsd(ip-8)+spsd(ip-9))
      psdsum=psdsum+amppsd*spsd(ip)/andim
      if(nTestOut.eq.6) then
         write(9,14009) 'psdsum',psdsum
         write(9,14009) 'spsd',spsd(ip)
      end if
 4998 continue
      if(nTestOut.eq.6) then
         close(9)
         call exit(0)
      end if 
      !write(6,9992)nend,ip0,dstart,zsize,zshock,clos,slos,
      !     ,   gammdd,delobs,dtfact,psdsig,psdsum
c      do 4999 ip=1,ndim
c      write(8,9345)ip,spsd(ip)
c 4999 continue
 9345 format(i10,1pe12.3)
c      ip0=randproto(0)*5000
      ip0=100
      it=0
 1000 continue
c     Set parameters of each cell at initial time
c     There are icells columns of cells, with jcells cells per column
c     The central cell is a Mach disk with low velocity; its observed radiation
c      is ignored, but it is an  important source of IC seed photons
      write(4,6667)
      write(5,6668)
      i=1
      j=0
      jzero=1
      do 999 nnn=1,(nend-1)
      do 998 jcnt=1,(6*nnn)
      j=j+1
      annn=nnn
      angc=60.0/(rad*annn)
      rcell(j)=2.0*rsize*annn
      xcell(j)=rcell(j)*cos(angc*(j-jzero))
      ycell(j)=rcell(j)*sin(angc*(j-jzero))
      cosph(j)=1.0
      sinph(j)=0.0
      zcol=rcell(j)/tanz
      imax(j)=2.0*zcol/zsize+0.01
      if(nTestOut.eq.2) write(9, 9219) j, imax(j)
      if(imax(j).lt.2)imax(j)=2
c     nouter(j) = approx. no. of cells between cell of interest and observer
      nouter(j)=imax(j)
      cosph(j)=xcell(j)/rcell(j)
      sinph(j)=ycell(j)/rcell(j)
      tanpsi(j)=rcell(j)/(zsvtex-zcol)
      psi(j)=datan(tanpsi(j))
      cospsi(j)=dcos(psi(j))
      sinpsi(j)=tanpsi(j)*cospsi(j)
  998 continue
      jzero=j+1
  999 continue
      
      if(nTestOut.eq.2) then
         close(9)
         call exit(0)
      end if

      xcell(jcells)=0.0
      ycell(jcells)=0.0
      rcell(jcells)=0.0
      cosph(jcells)=1.0
      sinph(jcells)=0.0
      zrf=zshock
      xrf=0.0
c
c     *** Set up Mach disk emission for earlier times ***
c
c     Compute time delay (no. of time steps) between when plasma passes Mach disk
c       and when it passes reference point of conical shock, in plasma frame
      zmd=zrf-zshock
      idelmd=zmd/(zsize/betd)
c      write(5,9595)mdmax,dstart,idelmd,ip0,zshock,rshock,zrf,zmd,
c     , zsize
c 9595 format('mdmax = ',i7,' dstart = ',i7,' idelmd = ',i7,
c     ,  '  ip0 = ',i7/'zshock, rshock, zrf, zmd, zsize ',
c     ,     1p5e10.2)
c     Velocity parameters of plasma in Mach disk
      betadx(1,jcells)=0.0
      betady(1,jcells)=0.0
      betadz(1,jcells)=1.0d0/3.0d0
      betamd=betadz(1,jcells)
      gammd=1.0d0/dsqrt(1.0d0-betamd**2)
      betad(1,jcells)=betamd
      gammad(1,jcells)=gammd
      dopref=gammdd/gammd
      i=1
c     Initialize some parameters for selecting magnetic field vector
      do 6129 j=1,jcells
      ididb(j)=0
      phi1(j)=0.0
 6129 bfrac(j)=0.0
c     Determine B vector of MD assuming random magnetic field orientation
c     Randomly select magnetic field direction for every 10th cell, then interpolate
c       inside loop to get direction for intermediate cells
c     First need to set up field direction for the more downstream cell
      phi2(jcells)=2.0*pi*randproto(0)
      costh=2.0*(randproto(0)-0.5)
      if(costh.ge.1.0)costh=0.999999
      if(costh.le.-1.0)costh=-0.999999
      theta2(jcells)=acos(costh)
c     Loop over mdmax cells that pass through Mach disk
      do 130 md=1,(mdmax-1)
      if(ididb(jcells).eq.1)go to 6130
      phi1(jcells)=2.0*pi*randproto(0)
      costh=2.0*(randproto(0)-0.5)
      if(costh.ge.1.0)costh=0.999999
      if(costh.le.-1.0)costh=-0.999999
      theta1(jcells)=acos(costh)
      sint1=sin(theta1(jcells))
      sint2=sin(theta2(jcells))
      bu1x(jcells)=sint1*cos(phi1(jcells))
      bu1y(jcells)=sint1*sin(phi1(jcells))
      bu1z(jcells)=costh
      bu2x(jcells)=sint2*cos(phi2(jcells))
      bu2y(jcells)=sint2*sin(phi2(jcells))
      bu2z(jcells)=cos(theta2(jcells))
      cosrot(jcells)=bu1x(jcells)*bu2x(jcells)+bu1y(jcells)*
     ,  bu2y(jcells)+bu1z(jcells)*bu2z(jcells)
      if(cosrot(jcells).ge.1.0)cosrot(jcells)=0.999999
      if(cosrot(jcells).le.-1.0)cosrot(jcells)=-0.999999
      angrot(jcells)=acos(cosrot(jcells))
      xsign=randproto(0)-0.5
      xsign=xsign/abs(xsign)
      if(xsign.lt.0.0)angrot(jcells)=angrot(jcells)-2.0*pi
      cu1x(jcells)=bu1y(jcells)*bu2z(jcells)-bu1z(jcells)*bu2y(jcells)
      cu1y(jcells)=-bu1x(jcells)*bu2z(jcells)+bu1z(jcells)*bu2x(jcells)
      cu1z(jcells)=bu1x(jcells)*bu2y(jcells)-bu1y(jcells)*bu2x(jcells)
      ididb(jcells)=1
 6130 call vecrot(bu1x(jcells),bu1y(jcells),bu1z(jcells),
     ,  cu1x(jcells),cu1y(jcells),cu1z(jcells),
     ,  (bfrac(jcells)*angrot(jcells)),bux,buy,buz)
      bfrac(jcells)=bfrac(jcells)+0.1
      if(bfrac(jcells).le.1.0)go to 6131
      theta2(jcells)=theta1(jcells)
      phi2(jcells)=phi1(jcells)
      bfrac(jcells)=0.0
      ididb(jcells)=0
 6131 continue
c     Compute B field components downstream of Mach disk shock
      idelay=dstart-mdmd+md-idelmd
      if(idelay.lt.1.or.idelay.gt.(ndim-ip0))
     ,  write(6,9299)idelay,i,jcells,md,idelmd,ip0,mdmax,dstart,
     ,  zshock,delobs,gamd,betd,clos
c      bavg=bave*sqrt(spsd(idelay+ip0))
      bavg=bave
      n0mean=n0ave*spsd(idelay+ip0)
      j=jcells
      n0(i,j)=etac*n0mean
      bx(i,j)=bavg*bux*etac
      by(i,j)=bavg*buy*etac
      bz(i,j)=bavg*buz
      bfield(j)=sqrt(bx(i,j)**2+by(i,j)**2+bz(i,j)**2)
      bfld=bfield(j)
      bmdx(md)=bx(i,j)
      bmdy(md)=by(i,j)
      bmdz(md)=bz(i,j)
      bmdtot(md)=bfld
c      write(6,9333)md,idelay,dstart,idelmd,bfld,bx(i,j),
c     ,  by(i,j),bz(i,j),bavg,n0mean,etac
c     Calculate the initial maximum electron energy in the Mach disk
      ball2=bux*bux+buy*buy+buz*buz
c     Next line relates this to direction of upstream B field relative to shock
c      gmax0(i,j)=gmaxmx*buz*buz/ball2
c      if(gmax0(i,j).lt.gmaxmn)gmax0(i,j)=gmaxmn
c     Next 3 lines assume a power-law distribution of gmax0, unrelated
c         to direction of B field
c     See http://mathworld.wolfram.com/RandomNumber.html for choosing
c      random numbers from a power-law distribution
c      xrand=rand(0)
c      if(pexp.eq.0.0)gmax0(i,j)=gmaxmx
c      if(pexp.lt.0.0)gmax0(i,j)=((gmaxmx**pexp-gmaxmn**pexp)*xrand+
c     ,  gmaxmn**pexp)**(1.0/pexp)
      gmax0(i,j)=gmaxmn
      gminmd=gmin
c      gminmd=1800.0
c      gmax0(i,j)=1.4*gminmd
c     Calculate energy distribution in the Mach disk cell
c     Compute energy density of photons in Mach disk, time delayed by 1 step
c       Ignore photons from dust torus since beaming is small in Mach disk
      mdi=md-1
      if(mdi.le.0)mdi=1
c     Value of SSC photon density for fast cooling case, from Sari and Esen (2001)
      uphmd=bmdtot(mdi)**2/(8.0*pi)*(sqrt(1.0+4.0*uratio)-1.0)/2.0
      ustob=8.0*pi*uphmd/(bfield(j))**2
      tlfact=cc2*bfield(j)**2*(1.0+ustob)
c     Effective length of Mach disk zsmd = loss time of highest-E
c       electrons (to energy glow where radiation is negligible) times flow velocity
      glow=10.0
      delt=(gmax0(i,j)-glow)/(tlfact*glow*gmax0(i,j))
      zsmd(md)=delt/(yr*3.26)*gammd*betamd
c     Use 0.99 instead of 1 to avoid singularity of N(gamma) at gmax0
      gmrat=0.99*gmax0(i,j)/gminmd
      gmratl=dlog10(gmrat)/32.0
      gmratm=dlog10(gminmd/glow)/11.0
c      write(6,9996)bfield(j),n0(i,j),gmax0(i,j),uphmd,ustob,
c     ,  tlfact,bmdtot(mdi),delt,gamb,glow
      do 124 ie=1,44
      if(ie.ge.12)egam(i,j,ie)=gminmd*10.0**(gmratl*(ie-12))
      if(ie.lt.12)egam(i,j,ie)=glow*10.0**(gmratm*(ie-1))
  123 ggam(ie)=egam(i,j,ie)
      t2=(gmax0(i,j)-ggam(ie))/(tlfact*ggam(ie)*gmax0(i,j))
      enofe(i,j,ie)=0.0
      eterm1=ggam(ie)/gminmd
      eterm2=1.0d0-ggam(ie)*t2*tlfact
      if(eterm2.lt.0.0)go to 5123
      if(ie.ge.12)enofe(i,j,ie)=n0(i,j)/((sen-1.0)*tlfact*
     ,  ggam(ie)**(sen+1.0))*(1.0d0-eterm2**(sen-1.0))
      if(ie.lt.12)enofe(i,j,ie)=n0(i,j)/((sen-1.0)*tlfact*
     ,  ggam(ie)**(sen+1.0))*(eterm1**(sen-1.0)-
     ,  eterm2**(sen-1.0))
c     Divide by delt since integral is over time
      delt=zsize*yr*3.26/(gammad(i,j)*betad(i,j))
      enofe(i,j,ie)=enofe(i,j,ie)/delt
      if(enofe(i,j,ie).lt.0.0)enofe(i,j,ie)=0.0
 5123 edist(ie)=enofe(i,j,ie)
      if((nTestOut.eq.7).and.(md>110000)) then
         delt_local = 0.0
         if(ie.eq.1) delt_local=delt
         write(9,14008) '5123',i,j,ie,md,delt_local,tlfact,n0(i,j),n0mean,etac,edist(ie)
      end if
c      aled=0.0
c      if(ie.gt.1.and.edist(ie).gt.0.0)aled=
c     ,  alog10(edist(ie)/edist(ie-1))/dlog10(ggam(ie)/ggam(ie-1))
c      write(6,9994)md,ie,bfield(j),uphmd,ggam(ie),edist(ie),aled
  124 continue
      if(nTestOut.eq.7) then
         cycle
      end if
      bperpp=bfld
      nuhi=1
      do 125 inu=1,40
      alfmdc(inu,md)=10.0
      fsscmd(inu,md)=0.0
      snu(inu)=nu(inu)
      restnu=nu(inu)
c     Synchrotron mean intensity for SSC calculation inside Mach disk
      ssabs=1.02e4*(sen+2.0)*akapnu(restnu)*bperpp/nu(inu)**2
      sstau=ssabs*parsec*rsize*svmd
      fsync(inu)=ajnu(restnu)*bperpp*rsize*svmd*(emfact*parsec)
      !if(nTestOut.eq.3) write(9, 9217) md, inu, fsync(inu)
      if(fsync(inu).gt.0.0)nuhi=inu
      if(sstau.lt.0.01)go to 126
      srcfn=fsync(inu)/sstau
      if(sstau.gt.5.0)fsync(inu)=srcfn
      if(sstau.gt.0.01.and.sstau.le.5.0)fsync(inu)=
     ,  srcfn*(1.0-exp(-sstau))
  126 continue
      !if(nTestOut.eq.3) write(9, 9217) md, inu, fsync(inu)
c     Now estimate synchrotron emission seen by other cells
c     Need to add a lower frequency to get spectral index of nu(1)
      if(inu.gt.1)go to 127
      fq1=0.98*nu(1)
      fsyn1=ajnu(fq1/dopref)*dopref**2*bperpp/dtfact
      absrb1=1.02e4*(sen+2.0)*akapnu(fq1/dopref)*bperpp/(fq1/dopref)**2
  127 continue
      fsynmd(inu,md)=ajnu(restnu/dopref)*dopref**2*bperpp/dtfact
      absorb(inu,md)=1.02e4*(sen+2.0)*akapnu(restnu/dopref)*
     ,   bperpp/nu(inu)**2
      alphmd(inu,md)=10.0
      abexmd(inu,md)=1.7
      if(fsynmd(inu,md).gt.0.0.and.fsyn1.gt.0.0)
     ,alphmd(inu,md)=-alog10(fsynmd(inu,md)/fsyn1)/alog10(restnu/fq1)
      if(absorb(inu,md).gt.0.0.and.absrb1.gt.0.0)
     ,abexmd(inu,md)=-alog10(absorb(inu,md)/absrb1)/alog10(restnu/fq1)
ccccccccc
c      write(6,9994)md,inu,restnu,fsync(inu),(restnu/dopref),
c     ,   fsynmd(inu,md),sstau,
c     ,   absorb(inu,md),alphmd(inu,md),bperpp,dopref
      fq1=restnu
      fsyn1=fsynmd(inu,md)
      absrb1=absorb(inu,md)
      ssseed(inu)=fsync(inu)
      !if(nTestOut.eq.3) write(9, 9216) md, inu, fsynmd(inu,md)
  125 continue
      do 128 inu=41,68
      snu(inu)=nu(inu)
      alphmd(inu,md)=10.0
      ssseed(inu)=0.0
      fsync(inu)=0.0
      fsynmd(inu,md)=0.0
      absorb(inu,md)=0.0
  128 continue
c     Calculate the SSC emission from the Mach disk for reference
c       relative Doppler factor dopref
      fq1=0.98*nu(7)
      betd=betamd
      gamd=gammd
      fssc1=ssc(fq1/dopref)*dopref**2/dtfact
      do 129 inu=7,68
      restnu=nu(inu)
      fsscmd(inu,md)=ssc(restnu/dopref)*dopref**2/dtfact
      alfmdc(inu,md)=10.0
      if(fsscmd(inu,md).gt.0.0.and.fssc1.gt.0.0)
     ,alfmdc(inu,md)=-alog10(fsscmd(inu,md)/fssc1)/alog10(restnu/fq1)
      fq1=restnu
      fssc1=fsscmd(inu,md)
      if(nTestOut.eq.3) write(9, 9215) md, inu, fsscmd(inu,md)
c      write(6,9911)md,j,dopref,bperpp,n0(i,j),
c     ,  restnu,ssseed(inu),fsynmd(inu,md),fsscmd(inu,md)
  129 continue
  130 continue
c
c     *** End Mach disk set-up ***
c
      if((nTestOut.eq.3).or.(nTestOut.eq.7)) then
         close(9)
         call exit(0)
      end if

c     After initial set-up, time loop starts here
    9 i=1
      ncells=0
      it=it+1
      md=mdmax
c
c     Start loop over all cells in first layer to set up physical parameters
c
      do 80 j=1,jcells
c
      do 81 inu=1,68
      fsynmd(inu,md)=0.0
      fsscmd(inu,md)=0.0
      fmdall(inu)=0.0
   81 continue
      ididg(j)=0
      zcell(i,j)=zshock-(rcell(j)-rsize)/tanz
      idelay=0.5+dstart+it-1+((zrf-zcell(i,j))
     ,   +(xrf-xcell(j))*betadd*slos/(1.0d0-betadd*clos))/zsize
      if(idelay.le.0.or.idelay.gt.(ndim-ip0))
     ,  write(5,9222)idelay,dstart,j,md,ip0,zshock,zrf,zcell(i,j),
     ,  xrf,xcell(j),betadd,zsize
c      write(6,9333)it,i,j,idelay,zcell(i,j),zref,xcell(j),betadd,zsize
c      bavg=bave*sqrt(spsd(idelay+ip0))
      bavg=bave
      n0mean=n0ave*spsd(idelay+ip0)
c     xxxxx
c      write(6,9333)it,i,j,idelay,zcell(i,j),rcell(j),xcell(j),
c     ,  ycell(j),spsd(idelay+ip0),bavg,n0mean
      phcell(j)=datan2(sinph(j),cosph(j))
c     Velocity vector of laminar component of pre-shock flow
      betupx=betaup*cosph(j)*sinpsi(j)
      betupy=betaup*sinph(j)*sinpsi(j)
      betupz=betaup*cospsi(j)
c     Velocity vector of the turbulent component of pre-shocked plasma
      phit=2.0*pi*randproto(0)
      costht=2.0d0*(randproto(0)-0.5d0)
      if(costht.gt.1.0d0)costht=0.99999
      if(costht.lt.-1.0d0)costht=-0.99999
      thetat=dacos(costht)
      sintht=dsin(thetat)
      betatx=betat*dcos(phit)*sintht
      betaty=betat*dsin(phit)*sintht
      betatz=betat*costht
      dotprd=betupx*betatx+betupy*betaty+betupz*betatz
      btparx=dotprd*betupx/betaup**2
      btpary=dotprd*betupy/betaup**2
      btparz=dotprd*betupz/betaup**2
      btprpx=betatx-btparx
      btprpy=betaty-btpary
      btprpz=betatz-btparz
c     Velocity vector of the pre-shock plasma including turbulent component
      betaux(j)=(betupx+btparx+btprpx/gamup)/(1.0d0+dotprd)
      betauy(j)=(betupy+btpary+btprpy/gamup)/(1.0d0+dotprd)
      betauz(j)=(betupz+btparz+btprpz/gamup)/(1.0d0+dotprd)
      betau(j)=dsqrt(betaux(j)**2+betauy(j)**2+betauz(j)**2)
      gammau(j)=1.0d0/dsqrt(1.0d0-betau(j)**2)
c     Unit vector of shock front at current position
      sx=-sinz*cosph(j)
      sy=-sinz*sinph(j)
      sz=cosz
c     Velocity vector downstream of shock + compression ratio of shock
      call vdcalc(betaux(j),betauy(j),betauz(j),sx,sy,sz,betadx(i,j),
     ,   betady(i,j),betadz(i,j),betad(i,j),gammad(i,j),eta)
cccccccccc
c      write(6,9994)i,j,betaux(j),betauy(j),betauz(j),betau(j),
c     ,  gammau(j),betadx(i,j),betady(i,j),betadz(i,j),betad(i,j),
c     ,  gammad(i,j),eta
c     Determine B vector of cell assuming random magnetic field orientation
   11 continue
c     For Mach disk, continue from previous calculation of B direction
      if(j.eq.jcells.and.ididb(j).eq.1)go to 6030
      if(j.eq.jcells)go to 6028
      if(ididb(j).eq.1)go to 6030
c     Randomly select magnetic field direction for every 10th cell, then interpolate
c       inside loop to get direction for intermediate cells
c     First need to set up field direction for the more downstream cell
      phi2(j)=2.0*pi*randproto(0)
      costh=2.0*(randproto(0)-0.5)
      if(costh.ge.1.0)costh=0.999999
      if(costh.le.-1.0)costh=-0.999999
      theta2(j)=acos(costh)
      sint1=sin(theta1(j))
      sint2=sin(theta2(j))
      bu1x(j)=sint1*cos(phi1(j))
      bu1y(j)=sint1*sin(phi1(j))
      bu1z(j)=cos(theta1(j))
      bu2x(j)=sint2*cos(phi2(j))
      bu2y(j)=sint2*sin(phi2(j))
      bu2z(j)=costh
      cosrot(j)=bu1x(j)*bu2x(j)+bu1y(j)*bu2y(j)+bu1z(j)*bu2z(j)
      if(cosrot(j).ge.1.0)cosrot(j)=0.999999
      if(cosrot(j).le.-1.0)cosrot(j)=-0.999999
      angrot(j)=acos(cosrot(j))
      xsign=randproto(0)-0.5
      xsign=xsign/abs(xsign)
      if(xsign.lt.0.0)angrot(j)=angrot(j)-2.0*pi
      cu1x(j)=bu1y(j)*bu2z(j)-bu1z(j)*bu2y(j)
      cu1y(j)=-bu1x(j)*bu2z(j)+bu1z(j)*bu2x(j)
      cu1z(j)=bu1x(j)*bu2y(j)-bu1y(j)*bu2x(j)
      ididb(j)=1
      if(phi1(j).ne.0.0)go to 6030
 6028 continue
c     Now set up field direction for cell just crossing the shock
      phi1(j)=2.0*pi*randproto(0)
      costh=2.0*(randproto(0)-0.5)
      if(costh.ge.1.0)costh=0.999999
      if(costh.le.-1.0)costh=-0.999999
      theta1(j)=acos(costh)
      sint1=sin(theta1(j))
      sint2=sin(theta2(j))
      bu1x(j)=sint1*cos(phi1(j))
      bu1y(j)=sint1*sin(phi1(j))
      bu1z(j)=costh
      bu2x(j)=sint2*cos(phi2(j))
      bu2y(j)=sint2*sin(phi2(j))
      bu2z(j)=cos(theta2(j))
      cosrot(j)=bu1x(j)*bu2x(j)+bu1y(j)*bu2y(j)+bu1z(j)*bu2z(j)
      if(cosrot(j).ge.1.0)cosrot(j)=0.999999
      if(cosrot(j).le.-1.0)cosrot(j)=-0.999999
      angrot(j)=acos(cosrot(j))
      xsign=randproto(0)-0.5
      xsign=xsign/abs(xsign)
      if(xsign.lt.0.0)angrot(j)=angrot(j)-2.0*pi
      cu1x(j)=bu1y(j)*bu2z(j)-bu1z(j)*bu2y(j)
      cu1y(j)=-bu1x(j)*bu2z(j)+bu1z(j)*bu2x(j)
      cu1z(j)=bu1x(j)*bu2y(j)-bu1y(j)*bu2x(j)
 6030 call vecrot(bu1x(j),bu1y(j),bu1z(j),cu1x(j),cu1y(j),cu1z(j),
     ,  (bfrac(j)*angrot(j)),bux,buy,buz)
c      buang1=atan2(buy,bux)
c      buang2=acos(buz)
      bfrac(j)=bfrac(j)+0.1
      if(bfrac(j).le.1.0)go to 6031
      theta1(j)=theta2(j)
      phi1(j)=phi2(j)
      bfrac(j)=0.0
      ididb(j)=0
 6031 continue
c     
c     Compute B field components downstream of shock in the plasma frame
c       by transforming the shock normal to the upstream plasma
c       frame, and compressing the component of B parallel to the shock
c
      bux=bavg*bux
      buy=bavg*buy
      buz=bavg*buz
c     Unit vector of shock normal at current position
      anx=cosz*cosph(j)
      any=cosz*sinph(j)
      anz=sinz
c     If Mach disk, compute separately
      if(j.eq.jcells)go to 12
c     Calculate upstream B field components parallel + perpendicular to shock front
      call bdcalc(betaux(j),betauy(j),betauz(j),betau(j),gammau(j),
     ,  anx,any,anz,bux,buy,buz,bparx,bpary,bparz,bprpx,bprpy,bprpz,
     ,  bpar,bprp)
      bx(i,j)=eta*bparx+bprpx
      by(i,j)=eta*bpary+bprpy
      bz(i,j)=eta*bparz+bprpz
c      write(6,9996)betaux(j),betauy(j),betauz(j),
c     ,  sx,sy,sz,bux,buy,buz,eta,
c     ,  bx(i,j),by(i,j),bz(i,j),bpar,bprp
c        write(6,9996)xcell(j),ycell(j),bx(i,j),by(i,j),bz(i,j),
c     ,    eta,bux,buy,buz,bparx,bpary,bparz,bpar,bprp
      go to 13
c     Set field of plasma in central cell, which is a Mach disk
   12 n0(i,j)=etac*n0mean
      bx(i,j)=bux*etac
      by(i,j)=buy*etac
      bz(i,j)=buz
      bfield(j)=sqrt(bx(i,j)**2+by(i,j)**2+bz(i,j)**2)
      bfld=bfield(j)
      bmdx(md)=bx(i,j)
      bmdy(md)=by(i,j)
      bmdz(md)=bz(i,j)
      bmdtot(md)=bfld
c     Calculate the initial maximum electron energy in the Mach disk
      ball2=bux*bux+buy*buy+buz*buz
c     Next line relates this to direction of upstream B field relative to shock
c      gmax0(i,j)=gmaxmx*buz*buz/ball2
c      if(gmax0(i,j).lt.gmaxmn)gmax0(i,j)=gmaxmn
c      gmax0(i,j)=gmaxmn
c     Next 3 lines assume a power-law distribution of gmax0, unrelated
c         to direction of B field
c     See http://mathworld.wolfram.com/RandomNumber.html for choosing
c      random numbers from a power-law distribution
c      xrand=rand(0)
c      if(pexp.eq.0.0)gmax0(i,j)=gmaxmx
c      if(pexp.lt.0.0)gmax0(i,j)=((gmaxmx**pexp-gmaxmn**pexp)*xrand+
c     ,  gmaxmn**pexp)**(1.0/pexp)
      gmax0(i,j)=gmaxmn
      gminmd=gmin
c      gminmd=1800.0
c      gmax0(i,j)=1.4*gminmd
c     Calculate energy distribution in the Mach disk cell
c     Compute energy density of photons in Mach disk, time delayed by 1 step
c       Ignore photons from dust torus since beaming is small in Mach disk
      mdi=md-1
      if(mdi.le.0)mdi=1
      uphmd=bmdtot(mdi)**2/(8.0*pi)*(sqrt(1.0+4.0*uratio)-1.0)/2.0
      ustob=8.0*pi*uphmd/(bfield(j))**2
      tlfact=cc2*bfield(j)**2*(1.0+ustob)
c     Effective length of Mach disk zsmd = loss time of highest-E
c       electrons (to energy glow where radiation is negligible) times flow velocity
      glow=10.0
      delt=(gmax0(i,j)-glow)/(tlfact*glow*gmax0(i,j))
      zsmd(md)=delt/(yr*3.26)*gammd*betamd
c      zsmd(md)=0.2*zsize
c     Use 0.99 instead of 1 to avoid singularity of N(gamma) at gmax0
      gmrat=0.99*gmax0(i,j)/gminmd
      gmratl=dlog10(gmrat)/32.0
      gmratm=dlog10(gminmd/glow)/11.0
c      write(6,9996)bfield(j),n0(i,j),gmax0(i,j),uphmd,ustob,
c     ,  tlfact,bmdtot(mdi),delt,gamb,glow
      do 1124 ie=1,44
      if(ie.ge.12)egam(i,j,ie)=gminmd*10.0**(gmratl*(ie-12))
      if(ie.lt.12)egam(i,j,ie)=glow*10.0**(gmratm*(ie-1))
 1123 ggam(ie)=egam(i,j,ie)
      t2=(gmax0(i,j)-ggam(ie))/(tlfact*ggam(ie)*gmax0(i,j))
      enofe(i,j,ie)=0.0
      eterm1=ggam(ie)/gminmd
      eterm2=1.0d0-ggam(ie)*t2*tlfact
      if(eterm2.lt.0.0)go to 5124
      if(ie.ge.12)enofe(i,j,ie)=n0(i,j)/((sen-1.0)*tlfact*
     ,  ggam(ie)**(sen+1.0))*(1.0d0-eterm2**(sen-1.0))
      if(ie.lt.12)enofe(i,j,ie)=n0(i,j)/((sen-1.0)*tlfact*
     ,  ggam(ie)**(sen+1.0))*(eterm1**(sen-1.0)-
     ,  eterm2**(sen-1.0))
c     Divide by delt since integral is over time
      delt=zsize*yr*3.26/(gammad(i,j)*betad(i,j))
      enofe(i,j,ie)=enofe(i,j,ie)/delt
      if(enofe(i,j,ie).lt.0.0)enofe(i,j,ie)=0.0
 5124 edist(ie)=enofe(i,j,ie)
      if((nTestOut.eq.7).and.(it.eq.1)) then
         delt_local = 0.0
         if(ie.eq.1) delt_local=delt
         write(9,14008) '5124',i,j,ie,delt_local,tlfact,n0(i,j),n0mean,etac,edist(ie)
      end if
 1124 continue
      bperpp=bfld
      nuhi=1
      do 1125 inu=1,40
      alfmdc(inu,md)=10.0
      fsscmd(inu,md)=0.0
      snu(inu)=nu(inu)
      restnu=nu(inu)
c     Synchrotron mean intensity for SSC calculation inside Mach disk
      ssabs=1.02e4*(sen+2.0)*akapnu(restnu)*bperpp/nu(inu)**2
      sstau=ssabs*parsec*rsize*svmd
c     fsync(inu) is the optically thin intensity
      fsync(inu)=ajnu(restnu)*bperpp*rsize*svmd*(emfact*parsec)
      if(fsync(inu).gt.0.0)nuhi=inu
      if(sstau.lt.0.01)go to 1126
c     Optically thick source function = emis. coef./abs. coef.
c     Note that the length of the path through the source cancels out
      srcfn=fsync(inu)/sstau
      if(sstau.gt.5.0)fsync(inu)=srcfn
      if(sstau.gt.0.01.and.sstau.le.5.0)fsync(inu)=
     ,  srcfn*(1.0-exp(-sstau))
 1126 continue
c     Now estimate synchrotron emission seen by other cells
c     Need to add a lower frequency to get spectral index of nu(1)
      if(inu.gt.1)go to 1127
      fq1=0.98*nu(1)
      fsyn1=ajnu(fq1/dopref)*dopref**2*bperpp/dtfact
      absrb1=1.02e4*(sen+2.0)*akapnu(fq1/dopref)*bperpp/(fq1/dopref)**2
 1127 continue
! TODO: In the next line, I matched the pattern of the previous dtfact occurrences instead
! of what the new temz.f does. Have sent an email to Marscher to confirm what he intended (Oct 2012)
      fsynmd(inu,md)=ajnu(restnu/dopref)*dopref**2*bperpp/dtfact
      absorb(inu,md)=1.02e4*(sen+2.0)*akapnu(restnu/dopref)*
     ,   bperpp/nu(inu)**2
      alphmd(inu,md)=10.0
      abexmd(inu,md)=1.7
      if(fsynmd(inu,md).gt.0.0.and.fsyn1.gt.0.0)
     ,alphmd(inu,md)=-alog10(fsynmd(inu,md)/fsyn1)/alog10(restnu/fq1)
      if(absorb(inu,md).gt.0.0.and.absrb1.gt.0.0)
     ,abexmd(inu,md)=-alog10(absorb(inu,md)/absrb1)/alog10(restnu/fq1)
c      write(5,9996)restnu,fsync(inu),fsynmd(inu,md),sstau,
c     ,   absorb(inu,md),alphmd(inu,md),bperpp
      fq1=restnu
      fsyn1=fsynmd(inu,md)
      absrb1=absorb(inu,md)
      ssseed(inu)=fsync(inu)
      if(nTestOut.eq.4) write(9,14001) 1125,i,j,md,inu,'fsync(inu)',fsync(inu)
      if(nTestOut.eq.4) write(9,14001) 1125,i,j,md,inu,'fsynmd(inu,md)',fsynmd(inu,md)
 1125 continue
      do 1128 inu=41,68
      snu(inu)=nu(inu)
      alphmd(inu,md)=10.0
      ssseed(inu)=0.0
      fsync(inu)=0.0
      fsynmd(inu,md)=0.0
      absorb(inu,md)=0.0
 1128 continue
c     Calculate the SSC emission from the Mach disk for reference
c       relative Doppler factor dopref
      fq1=0.98*nu(7)
      betd=betamd
      gamd=gammd
      fssc1=ssc(fq1/dopref)*dopref**2/dtfact
      do 1129 inu=7,68
      restnu=nu(inu)
      fsscmd(inu,md)=ssc(restnu/dopref)*dopref**2/dtfact
      alfmdc(inu,md)=10.0
      if(fsscmd(inu,md).gt.0.0.and.fssc1.gt.0.0)
     ,alfmdc(inu,md)=-alog10(fsscmd(inu,md)/fssc1)/alog10(restnu/fq1)
c      write(5,9996)restnu,fsscmd(inu,md),alfmdc(inu,md)
      fq1=restnu
      fssc1=fsscmd(inu,md)
 1129 continue
      go to 88
   13 continue
c     Calculate component of magnetic field that is perpendicular to
c       the aberrated line of sight in the plasma frame
c     Line-of-sight unit vector
      slx=slos
      sly=0.0
      slz=clos
      call bcalc(betadx(i,j),betady(i,j),betadz(i,j),betau(j),gammau(j),
     ,  slx,sly,slz,bx(i,j),by(i,j),bz(i,j),dum1,dum2,dum3,dum4,
     ,  dum5,dum6,dum7,bperp(j))
      bfield(j)=sqrt(bx(i,j)**2+by(i,j)**2+bz(i,j)**2)
      n0(i,j)=eta*n0mean
c     Calculate the initial maximum electron energy in each cell
c       from the ratio of B(perp. to shock) to B_total in shock frame
      gmax0(i,j)=gmaxmx*(bprp*bprp/(bpar*bpar+bprp*bprp))
      if(gmax0(i,j).lt.gmaxmn)gmax0(i,j)=gmaxmn
c     Next 3 lines assume a power-law distribution of gmax0, unrelated
c         to direction of B field
c     See http://mathworld.wolfram.com/RandomNumber.html for choosing
c      random numbers from a power-law distribution
c      xrand=randproto(0)
c      if(pexp.eq.0.0)gmax0(i,j)=gmaxmx
c      if(pexp.lt.0.0)gmax0(i,j)=((gmaxmx**pexp-gmaxmn**pexp)*xrand+
c     ,  gmaxmn**pexp)**(1.0/pexp)
      do 15 ig=1,43
   15 if(gmax0(i,j).le.gcnt(ig+1).and.gmax0(i,j).ge.gcnt(ig))
     , igcnt(ig)=igcnt(ig)+1
   80 continue
c
c     End initial loop over all cells in first layer
c
c     Start another loop over all cells in first layer to calculate the
c       electron energy distribution and the emission
c
   88 continue
      do 100 j=1,jcells-1
      ncells=ncells+1
      !if(nTestOut.eq.2) write(9,9220) ncells, j
      emisco=0.0
      ecflux=0.0
      zcell(i,j)=zshock-(rcell(j)-rsize)/tanz
c     Determine Doppler factor relative to the observer
      bdx=betadx(i,j)
      bdy=betady(i,j)
      bdz=betadz(i,j)
      betacs=bdx*slos+bdz*clos
      delta(i,j)=1.0d0/(gammad(i,j)*(1.0d0-betacs))
c      write(4,9994)i,j,bdx,bdy,bdz,betacs,gammad(i,j),delta(i,j)
      ecflux=0.0
      bperpp=bperp(j)
      bfld=bfield(j)
c     Determine velocity of the cell plasma relative to the MD plasma
      bmparx=betamd*bdx*bdz/betad(i,j)**2
      bmpary=betamd*bdy*bdz/betad(i,j)**2
      bmparz=betamd*bdz*bdz/betad(i,j)**2
      bmprpx=-bmparx
      bmprpy=-bmpary
      bmprpz=betamd-bmparz
      betamx(j)=(bdx-bmparx-bmprpx/gammad(i,j))/(1.0d0-betamd*bdz)
      betamy(j)=(bdy-bmpary-bmprpy/gammad(i,j))/(1.0d0-betamd*bdz)
      betamz(j)=(bdz-bmparz-bmprpz/gammad(i,j))/(1.0d0-betamd*bdz)
      betamr(j)=dsqrt(betamx(j)**2+betamy(j)**2+betamz(j)**2)
      gamamr(j)=1.0d0/dsqrt(1.0d0-betamr(j)**2)
      zcl=zcell(i,j)
      zrel=zshock-0.5*zsize-zcl
      dmd(j)=sqrt(zrel**2+rcell(j)**2)
c     Calculate angle between MD seed photon and scattered l.o.s. photon
      call scatcs(bdx,bdy,bdz,clos,0.0d0,slos,xcell(j),ycell(j),
     ,   -zrel,cscat)
c     Determine Doppler factor of the MD emission in the cell's frame
      dmdp=sqrt(rcell(j)**2+(zrel/gamamr(j))**2)
      cosmd=(betamx(j)*xcell(j)+betamy(j)*ycell(j)+
     ,  betamz(j)*zrel/gamamr(j))/(betamr(j)*dmdp)
      tanmd=dsqrt(1.0d0/cosmd**2-1.0d0)
      angm=datan(tanmd)
      tanmd=dtan(2.0*datan(dtan(0.5d0*angm)/
     ,  (gamamr(j)*(1.0d0+betamr(j)))))
      cosmd=1.0d0/dsqrt(1.0d0+tanmd**2)
      deltmd(j)=1.0d0/(gamamr(j)*(1.0d0-betamr(j)*cosmd))
c     Determine time step of Mach disk seed photons from light-travel delay
      delcor=deltmd(j)/dopref
      mdmid=mdmd-(((zrf-zcl)*clos+(xrf-
     ,   xcell(j))*slos)+dmd(j))/(dtfact*zsize)
      md1=mdmid-mdrang
      if(md1.lt.1)md1=1
      if(mdmid.lt.0.or.md1.gt.mdmax)write(6,9736)it,i,j,mdmid,
     , md1,idelay,zrf,zcl,xrf,xcell(j),dmd(j),dtfact,zsize,clos,slos
      md2=mdmid+mdrang
c      write(7,9191)md1,md2,mdmid
c 9191 format(3i10)
      if(md2.gt.mdmax)md2=mdmax
      if(md1.gt.md2)md1=md2
      amdrng=md2-md1+1.0
      do 3145 inu=1,68
 3145 fmdall(inu)=0.0
      do 2147 md=md1,md2
      do 3146 inu=1,68
      syseed(inu)=0.0
      scseed(inu)=0.0
      tauxmd(inu)=0.0
 3146 ssseed(inu)=0.0
c     Mach disk's B field in frame of cell plasma transverse to cell's
c       line of sight to Mach disk (for synchrotron calculation)
      sx=rcell(j)*cosph(j)
      sy=rcell(j)*sinph(j)
      sz=zrel
      call bcalc(betamx(j),betamy(j),betamz(j),betamr(j),gamamr(j),
     ,  sx,sy,sz,bmdx(md),bmdy(md),bmdz(md),dum1,dum2,dum3,dum4,
     ,  dum5,dum6,dum7,bmperp)
c     Apply as a correction factor to previous estimate of B_perpendicular
      bpcorr=bmperp/bmdtot(md)
c      write(5,9996)dmd(j),rcell(j),zrel,deltmd(j),gamamr(j),betamr(j),
c     , tanmd,cosmd,delcor,angm
c     , bmdx(md),bmdy(md),bmdz(md),betamx(j),
c     , betamy(j),betamz(j),bpcorr
c
c     Calculate emission from various processes
c
c     Calculate synchrotron flux from Mach disk in frame of cell plasma
c
      nuhi=1
      do 146 inu=1,40
      if(fsynmd(inu,md).eq.0.0.or.alphmd(inu,md).gt.9.0)
     ,   go to 145
c     Synchrotron mean intensity for inverse Compton calculation
      ssabs=absorb(inu,md)*bpcorr**(0.5*(sen+2.0))*
     ,  delcor**(abexmd(inu,md))
      path=2.0*rsize*svmd/(tanmd*cosmd)
      if(path.gt.zsmd(md))path=zsmd(md)
      ajofnu=fsynmd(inu,md)*bpcorr**(1.0+alphmd(inu,md))*
     ,  delcor**(2.0+alphmd(inu,md))
      ajofnu=ajofnu*zsmd(md)*(emfact*parsec)
      syseed(inu)=ajofnu
      srcfn=ajofnu/sstau
      if(sstau.gt.5.0)syseed(inu)=srcfn
      if(sstau.gt.0.1.and.sstau.le.5.0)syseed(inu)=
     ,  srcfn*(1.0-exp(-sstau))
      syseed(inu)=syseed(inu)*pi*rsize**2*vmd/(dmd(j)**2)
c     Now estimate exponential attenuation of MD seed photons by
c      synchrotron self-absorption in intervening cells
      tauexp=1.02e4*(sen+2.0)*akapnu(nu(inu))*
     ,  bperp(j)/nu(inu)**2*parsec*dmd(j)
      tauxmd(inu)=tauexp
      if(tauexp.gt.15.0)syseed(inu)=0.0
      if(tauexp.le.15.0)syseed(inu)=syseed(inu)/exp(tauexp)
  145 scseed(inu)=0.0
      if(inu.lt.16)fmdall(inu)=fmdall(inu)+syseed(inu)/
     ,  (amdrng/dtfact)
c      write(6,9996)snu(inu),restnu,ssabs,sstau,ajofnu,
c     ,  fsynmd(inu,md),bpcorr,delcor,alphmd(inu,md),syseed(inu)
  146 continue
c
c     Calculate inverse Compton flux from Mach disk in frame of cell plasma
c
      do 147 inu=16,68
c     Inverse Compton mean intensity from Mach disk for 2nd-order
c       inverse Compton calculation in cells
      scseed(inu)=fsscmd(inu,md)*
     ,  delcor**(2.0+alfmdc(inu,md))
      scseed(inu)=scseed(inu)*zsmd(md)*pi*rsize**2*
     ,  (amjy*parsec)*vmd/dmd(j)**2
      if(tauxmd(inu).gt.15.0)scseed(inu)=0.0
      if(tauxmd(inu).le.15.0)scseed(inu)=scseed(inu)/
     ,exp(tauxmd(inu))
      fmdall(inu)=fmdall(inu)+(syseed(inu)+scseed(inu))/
     ,  (amdrng/dtfact)
c      write(6,9996)snu(inu),restnu,tauxmd(inu),syseed(inu),
c     ,  scseed(inu)
  147 continue
 2147 continue
      nuhi=0
      do 2148 inu=1,68
      ssseed(inu)=fmdall(inu)
 2148 if(fmdall(inu).gt.0.0)nuhi=inu
c     Calculate seed photon energy density energy loss calculation
      usdmd=0.0
      do 149 inu=2,nuhi
c     Only calculate up to Klein-Nishina limit of gmin
      epslon=6.63e-27*snu(inu)/(gmin*emc2)
      if(epslon.gt.1.0)go to 150
      if(fmdall(inu).le.SMALL_FMDALL.or.fmdall(inu-1).le.SMALL_FMDALL)go to 148
      aaa=alog10(fmdall(inu-1)/fmdall(inu))/alog10(snu(inu)/snu(inu-1))
      usdmd=usdmd+0.5/c/(1.0-aaa)*fmdall(inu-1)*snu(inu-1)*
     ,  ((snu(inu)/snu(inu-1))**(1.0-aaa)-1.0)
  148 continue
  149 continue
  150 continue
c     Calculate energy distribution in the cell
      id=(zcell(i,j)-zshock)/zsize+nzid
      if(id.lt.1)id=1
      if(id.gt.(icells+nzid))write(6,9992)i,j,id,zcell(i,j),zshock,
     ,  zsize
      zdist=zdist0+(id-nend/zrat)*zsize+zshock
c     Calculate min & max angles of dust torus in plasma frame
      dphi1=asin(dtrad/sqrt(zdist**2+dtdist**2))
      dphi2=atan(dtdist/zdist)
      dth1=dphi2-dphi1
      dth2=dphi2+dphi1
      csth1=dcos(dth1)
      csth2=dcos(dth2)
      dcsth1=-(csth1-betad(i,j))/(1.0d0-csth1*betad(i,j))
      dcsth2=-(csth2-betad(i,j))/(1.0d0-csth2*betad(i,j))
      csang=0.5*(dcsth1+dcsth2)
c     Doppler factor of dust torus emission in frame of cells, used
c       to estimate frequency of peak intensity in plasma frame
      tdel=gammad(i,j)*(1.0d0-betad(i,j)*0.5d0*(csth1+csth2))
c     Correct energy density of seed photons from dust for Doppler factor
c       that include turbulent component of velocity; approximates that
c       motion is directed parallel to jet axis; also include first-order
c       dependence on scattering angle in Klein-Nishina cross-section
      useedr=useed(id)*(tdel/tdelr(id))**2/(1.0d0+betad(i,j)*csang)
      ustob=8.0*pi*(useedr+usdmd)/(bfield(j))**2
      tlfact=cc2*bfield(j)**2*(1.0+ustob)
      tlf1(j)=tlfact
      delt=zsize*yr*3.26/(gammad(i,j)*betad(i,j))
      gamb=gmax0(i,j)/(1.0+tlfact*gmax0(i,j)*delt)
      glow=gmin/(1.0+tlfact*gmin*delt)
      gmrat=0.99*gmax0(i,j)/gmin
      gmratl=dlog10(gmrat)/32.0
      gmratm=dlog10(gmin/glow)/11.0
      ibreak=0
      do 90 ie=1,44
      egam(i,j,ie)=gmin*10.0**(gmratl*(ie-12))
      if(ie.lt.12)egam(i,j,ie)=glow*10.0**(gmratm*(ie-1))
      if(ie.eq.1)egam(i,j,ie)=glow+0.2*(glow*10.0**(gmratm)-glow)
      if(ibreak.eq.1.or.egam(i,j,ie).lt.gamb)go to 89
      egam(i,j,ie)=gamb
      ibreak=1
   89 ggam(ie)=egam(i,j,ie)
      tloss=(gmax0(i,j)-ggam(ie))/(tlfact*ggam(ie)*gmax0(i,j))
      t2=dmin1(tloss,delt)
      enofe(i,j,ie)=0.0
      eterm1=ggam(ie)/gmin
      eterm2=1.0d0-ggam(ie)*t2*tlfact
      if(eterm2.lt.0.0)go to 5089
      if(ie.ge.12)enofe(i,j,ie)=n0(i,j)/((sen-1.0)*tlfact*
     ,  ggam(ie)**(sen+1.0))*(1.0d0-eterm2**(sen-1.0))
      if(ie.lt.12)enofe(i,j,ie)=n0(i,j)/((sen-1.0)*tlfact*
     ,  ggam(ie)**(sen+1.0))*(eterm1**(sen-1.0)-
     ,  eterm2**(sen-1.0))
c     Divide by delt since integral is over time
      enofe(i,j,ie)=enofe(i,j,ie)/delt
      if(enofe(i,j,ie).lt.0.0)enofe(i,j,ie)=0.0
 5089 edist(ie)=enofe(i,j,ie)
      if((nTestOut.eq.7).and.(it.eq.1)) then
         delt_local = 0.0
         if(ie.eq.1) delt_local=delt
         write(9,14008) '5089',i,j,ie,delt_local,tlfact,n0(i,j),n0mean,etac,edist(ie)
      end if
   90 continue
      delt=zsize*3.26*yr/(gammad(i,j)*betad(i,j))
      gammax(i,j)=gmax0(i,j)/(1.0d0+tlfact*delt*gmax0(i,j))
      gammin(i,j)=gmin/(1.0d0+tlfact*delt*gmin)
      emold=0.0
c     Start loop over observed frequencies
      ithin=0
      do 93 inu=1,22
      dustnu(inu)=dustf(id,inu)
      hnukt=4.8e-11*dustnu(inu)/tdust
      hnuktr=4.8e-11*dustnu(inu)*(tdel/tdelr(id))/tdust
      dusti(inu)=dustii(id,inu)
      if(hnuktr.gt.60.0)go to 5093
c   93 write(6,9994)id,inu,dustnu(inu),dustii(id,inu),dusti(inu),
c     ,    tdel,tdelr(id),tdust,hnukt,hnuktr
      if(dusti(inu).gt.1.0e-30)dusti(inu)=dusti(inu)*
     ,  (tdel/tdelr(id))**3*(1.0-exp(hnukt))/(1.0-exp(hnuktr))
      go to 5094
 5093 dusti(inu)=0.0
 5094 continue
   93 continue
      do 95 inu=1,68
      restnu=nu(inu)*zred1/delta(i,j)
      snu(inu)=nu(inu)
      specin=0.0001
      emisco=0.0
      fsync2(inu)=0.0
      ssabs=0.0
      ecflux=0.0
      sscflx=0.0
      flsync(i,j,inu)=0.0
      flcomp(i,j,inu)=0.0
      flec(i,j,inu)=0.0
      flssc(i,j,inu)=0.0
      flux(i,j,inu)=0.0
      if(inu.gt.44)go to 92
      emisco=ajnu(restnu)*bperp(j)*delta(i,j)**2
      fsync2(inu)=emisco*zsize*(emfact*parsec)
      fsnoab=fsync2(inu)
      if(inu.eq.1)go to 91
      if(emold.gt.0.0.and.fsync2(inu).gt.0.0)specin=
     ,  alog10(emold/fsync2(inu))/alog10(nu(inu)/nu(inu-1))
   91 continue
      if(ithin.eq.1)go to 92
      ssabs=1.02e4*(sen+2.0)*akapnu(restnu)*bperp(j)/
     ,   nu(inu)**2/delta(i,j)
      ssabs=ssabs*parsec*zsize*delta(i,j)
c     Attenuate Mach disk emission from synchrotron absorption on way to cell
      ssabsm=ssabs*dmd(j)/zsize
      if(ssabsm.ge.10.0)ssseed(inu)=0.0
      if(ssabsm.lt.10.0.and.ssabsm.gt.0.1)ssseed(inu)=
     ,     ssseed(inu)/exp(ssabsm)
c     Return to absorption within cell
c     Use rsize instead of zsize because of aberration
      ssabs=ssabs*rsize/zsize
      if(ssabs.le.(0.1/ancol))ithin=1
      srcfn=fsync2(inu)/ssabs
      if(ssabs.gt.5.0)fsync2(inu)=srcfn
      if(ssabs.gt.0.1.and.ssabs.le.5.0)fsync2(inu)=
     ,  srcfn*(1.0-exp(-ssabs))
      tauexp=nouter(j)*ssabs
      if(rcell(j).gt.(0.98*rbound).and.xcell(j).le.0.0)
     ,         tauexp=0.0
      if(thlos.eq.0.0)tauexp=nouter(j)*ssabs
c      if(tauexp.gt.15.0)fsync2(inu)=0.0
c      if(tauexp.le.15.0)fsync2(inu)=fsync2(inu)/exp(tauexp)
      specin=alpha
   92 continue
      flsync(i,j,inu)=fsync2(inu)*(volc/zsize)*zred1/
     ,  (1.0e18*amjy*dgpc**2)*fgeom
      flux(i,j,inu)=flsync(i,j,inu)
c xxxxxxx
c      if(inu.eq.20)write(6,9994)i,j,nu(inu),flsync(i,j,inu),specin,
c     ,  bperp(j),delta(i,j),ggam(43),edist(43)
      betd=betad(i,j)
      gamd=gammad(i,j)
      if(inu.eq.1)call polcalc(bfield(j),bx(i,j),by(i,j),bz(i,j),
     ,    clos,slos,chipol)
      if(specin.lt.alpha)specin=alpha
      poldeg=(specin+1.0)/(specin+5.0/3.0)
      if(ssabs.gt.1.0)poldeg=3.0/(12.0*specin+19)
      fpol(i,j,inu)=poldeg*flsync(i,j,inu)
      pq(i,j,inu)=fpol(i,j,inu)*cos(2.0*chipol)
      pu(i,j,inu)=fpol(i,j,inu)*sin(2.0*chipol)
c      if(i.eq.1.and.inu.eq.12)write(6,9992)it,i,j,bfield(j),
c     ,    bx(i,j),by(i,j),bz(i,j),clos,slos,pq(i,j,inu),
c     ,    pu(i,j,inu),fpol(i,j,inu),(chipol*rad)
c      if(inu.eq.20)write(6,9989)i,j,inu,nu(inu),restnu,n0(i,j),
c     ,  bperp(j),delta(i,j),ggam(40),ggam(43),edist(40),edist(43)
      if(restnu.lt.1.0e14)go to 94
      spxec=0.0001
      spxssc=0.0001
      betd=betad(1,jcells)
      gamd=gammad(1,jcells)
      ecdust_local = ecdust(restnu)
      ecflux=ecdust_local*3.086*volc*zred1*delta(i,j)**2/dgpc**2*
     ,   fgeom
      if((nTestOut.eq.1.).and.(it.eq.1)) write(9, 9218) 92, i, j, inu, ecdust_local, delta(i,j)
c      if(inu.eq.29)write(6,9973)i,j,ecflux,delta(i,j),dusti(1),
c     ,  dusti(22),ggam(1),edist(1),ggam(44),edist(44)
 9973 format('EC ',2i6,1p8e11.3)
      sscflx=ssc(restnu)*3.086*volc*zred1*delta(i,j)**2/dgpc**2*
     ,   fgeom
      taupp=0.0
      if(nu(inu).lt.1.0e22)go to 99
c     Pair production opacity calculation
c     Expression for anumin includes typical interaction angle
      anumin=(1.24e20/(nu(inu)*zred1))*1.24e20*(2.0*gammad(i,j))**2
      alnumn=alog10(anumin)
      inumin=(alnumn-10.0)*4+1
      if(inumin.gt.40)go to 99
      if(anumin.gt.nu(inumin))anumin=nu(inumin)
      do 97 iinu=inumin,40
      bp=sqrt(1.0-anumin/nu(iinu))
      if(bp.le.0.0)bp=0.001
      xsecp1=1.25e-25*(1.0-bp*bp)*((3.0-bp**4)*
     *  alog((1.0+bp)/(1.0-bp))-2.0*bp*(2.0-bp*bp))
      bp=sqrt(1.0-anumin/nu(iinu+1))
      if(bp.le.0.0)bp=0.001
      xsecp2=1.25e-25*(1.0-bp*bp)*((3.0-bp**4)*
     *  alog((1.0+bp)/(1.0-bp))-2.0*bp*(2.0-bp*bp))
      taupp=taupp+0.5*(phots(iinu)*xsecp1+phots(iinu+1)*xsecp2)*
     ,  (nu(iinu+1)-nu(iinu))*nouter(j)*zsize*parsec
   97 continue
c      write(5,9996)nu(inu),taupp,anumin,
c     ,   phots(inumin),nu(inumin+10),phots(inumin+10)
      if(taupp.lt.0.1)go to 99
      if(taupp.lt.10.0)go to 98
      ecflux=0.0
      sscflx=0.0
      go to 99
   98 ecflux=ecflux/exp(taupp)
      sscflx=sscflx/exp(taupp)
   99 flec(i,j,inu)=ecflux
      flssc(i,j,inu)=sscflx
      flcomp(i,j,inu)=ecflux+sscflx
      flux(i,j,inu)=flsync(i,j,inu)+flcomp(i,j,inu)
      if((nTestOut.eq.1.).and.(it.eq.1)) write(9, 9218) 99, i, j, inu, ecflux, delta(i,j)
      if(emeold.gt.0.0.and.ecflux.gt.0.0)spxec=
     ,  alog10(emeold/ecflux)/alog10(nu(inu)/nu(inu-1))
      if(emsold.gt.0.0.and.sscflx.gt.0.0)spxssc=
     ,  alog10(emsold/sscflx)/alog10(nu(inu)/nu(inu-1))
   94 continue
cc      if(inu.eq.53)write(5,9994)i,j,restnu,n0(i,j),gammax(i,j),
cc     ,   gammin(i,j),ecflux,sscflx,syseed(20),scseed(35)
      emold=fsnoab
      emeold=ecflux
      emsold=sscflx
c      if(nu(inu).eq.1.0e20)
c     ,write(5,9994)i,j,restnu,emisco,ecflux,flux(i,j,inu),bperpp,
c     ,  gammax(i-1,j),gammin(i-1,j),gammax(i,j),gammin(i,j)
c      if(inu.eq.21)write(5,9992)i,j,inu,ecflux,dusti(1),dusti(4),
c     ,  dusti(8),dusti(12),dusti(16),dusti(20),dusti(22)
   95 continue
ccc    Test
cc      if(it.eq.50.or.it.eq.140)write(5,9989)i,j,mdd(j),bperp(j),
cc     ,  n0(i,j),gammax(i,j),gammin(i,j),flec(i,j,52),flssc(i,j,52)
   96 continue
      icelmx(j)=1
  100 continue
c
c     End loop over first layer of cells
c
c      go to 2200
      if(icells.eq.1)go to 2200
      i=i+1
      istart=i
c     Start loop over downstream cells
      do 200 j=1,(jcells-1)
c     iend, to be computed, is the last slice of cells with energetic electrons
      iend=1
      ididg(j)=0
      do 200 i=istart,imax(j)
      icelmx(j)=i
      ncells=ncells+1
      !if(nTestOut.eq.2) write(9,9221) ncells, j, i 
      if(it.gt.1)go to 110
c
c     *** Initial set-up of downstream cells; skip after 1st time step
c
      zcell(i,j)=(i-1)*zsize+zshock-(rcell(j)-rsize)/tanz
c     Set up physical parameters of downstream cells at first time step
c     Cells farther downstream contain plasma ejected earlier
      idelay=0.5+dstart+it-1+((zrf-zcell(i,j))
     ,   +(xrf-xcell(j))*betadd*slos/(1.0d0-betadd*clos))/zsize
      if(idelay.lt.1.or.idelay.gt.(ndim-ip0))
     ,  write(5,9223)idelay,dstart,j,md,ip0,zshock,zrf,zcell(i,j),
     ,  xrf,xcell(j),betad(i,j),zsize
c      write(6,9333)it,i,j,idelay,zcell(i,j),zref,xcell(j),betadd,zsize
c      bavg=bave*sqrt(spsd(idelay+ip0))
      bavg=bave
      n0mean=n0ave*spsd(idelay+ip0)
c      write(6,9333)it,i,j,idelay,zcell(i,j),xcell(j),
c     ,  spsd(idelay+ip0),bavg,n0mean
c     Velocity vector of laminar component of pre-shock flow
      betupx=betaup*cosph(j)*sinpsi(j)
      betupy=betaup*sinph(j)*sinpsi(j)
      betupz=betaup*cospsi(j)
c     Velocity vector of the turbulent component of pre-shocked plasma
      phit=2.0*pi*randproto(0)
      costht=2.0d0*(randproto(0)-0.5d0)
      if(costht.gt.1.0d0)costht=0.99999
      if(costht.lt.-1.0d0)costht=-0.99999
      thetat=dacos(costht)
      sintht=dsin(thetat)
      betatx=betat*dcos(phit)*sintht
      betaty=betat*dsin(phit)*sintht
      betatz=betat*costht
      dotprd=betupx*betatx+betupy*betaty+betupz*betatz
      btparx=dotprd*betupx/betaup**2
      btpary=dotprd*betupy/betaup**2
      btparz=dotprd*betupz/betaup**2
      btprpx=betatx-btparx
      btprpy=betaty-btpary
      btprpz=betatz-btparz
c     Velocity vector of the pre-shock plasma including turbulent component
      betaux(j)=(betupx+btparx+btprpx/gamup)/(1.0d0+dotprd)
      betauy(j)=(betupy+btpary+btprpy/gamup)/(1.0d0+dotprd)
      betauz(j)=(betupz+btparz+btprpz/gamup)/(1.0d0+dotprd)
      betau(j)=dsqrt(betaux(j)**2+betauy(j)**2+betauz(j)**2)
      gammau(j)=1.0d0/dsqrt(1.0d0-betau(j)**2)
c     Unit vector of shock front at current position
      sx=-sinz*cosph(j)
      sy=-sinz*sinph(j)
      sz=-cosz
c     Velocity vector downstream of shock + compression ratio of shock
      call vdcalc(betaux(j),betauy(j),betauz(j),sx,sy,sz,betadx(i,j),
     ,   betady(i,j),betadz(i,j),betad(i,j),gammad(i,j),eta)
      betacs=betadx(i,j)*slos+betadz(i,j)*clos
      delta(i,j)=1.0d0/(gammad(i,j)*(1.0d0-betacs))
c      write(6,9994)i,j,betaux(j),betauy(j),betauz(j),betau(j),
c     ,  gammau(j),betadx(i,j),betady(i,j),betadz(i,j),betad(i,j),
c     ,  gammad(i,j),delta(i,j),eta
  103 continue
c     Randomly select magnetic field direction for every 10th cell, then interpolate
c       inside loop to get direction for intermediate cells
c     First need to initialize values to maintain continuity with it=1,i=1 values
      if(i.gt.2)go to 6229
      ididbd(j)=ididb(j)
      bfracd(j)=bfrac(j)
      thet1d(j)=theta1(j)
      thet2d(j)=theta2(j)
      phi1d(j)=phi1(j)
      phi2d(j)=phi2(j)
      angrtd(j)=angrot(j)
      cosrtd(j)=cos(angrtd(j))
      bu1xd(j)=bu1x(j)
      bu1yd(j)=bu1y(j)
      bu1zd(j)=bu1z(j)
      bu2xd(j)=bu2x(j)
      bu2yd(j)=bu2y(j)
      bu2zd(j)=bu2z(j)
      cu1xd(j)=cu1x(j)
      cu1yd(j)=cu1y(j)
      cu1zd(j)=cu1z(j)
 6229 if(ididbd(j).eq.1)go to 6230
      phi2d(j)=2.0*pi*randproto(0)
      costh=2.0*(randproto(0)-0.5)
      if(costh.ge.1.0)costh=0.999999
      if(costh.le.-1.0)costh=-0.999999
      thet2d(j)=acos(costh)
      if(i.le.2)go to 6228
      sint1=sin(thet1d(j))
      sint2=sin(thet2d(j))
      bu1xd(j)=sint1*cos(phi1d(j))
      bu1yd(j)=sint1*sin(phi1d(j))
      bu1zd(j)=cos(thet1d(j))
      bu2xd(j)=sint2*cos(phi2d(j))
      bu2yd(j)=sint2*sin(phi2d(j))
      bu2zd(j)=costh
      cosrtd(j)=bu1xd(j)*bu2xd(j)+bu1yd(j)*bu2yd(j)+bu1zd(j)*bu2zd(j)
      if(cosrtd(j).ge.1.0)cosrtd(j)=0.999999
      if(cosrtd(j).le.-1.0)cosrtd(j)=-0.999999
      angrtd(j)=acos(cosrtd(j))
      xsign=randproto(0)-0.5
      xsign=xsign/abs(xsign)
      if(xsign.lt.0.0)angrtd(j)=angrtd(j)-2.0*pi
      cu1xd(j)=bu1yd(j)*bu2zd(j)-bu1zd(j)*bu2yd(j)
      cu1yd(j)=-bu1xd(j)*bu2zd(j)+bu1zd(j)*bu2xd(j)
      cu1zd(j)=bu1xd(j)*bu2yd(j)-bu1yd(j)*bu2xd(j)
      ididbd(j)=1
      if(i.gt.2)go to 6230
 6228 phi1d(j)=2.0*pi*randproto(0)
      costh=2.0*(randproto(0)-0.5)
      if(costh.ge.1.0)costh=0.999999
      if(costh.le.-1.0)costh=-0.999999
      thet1d(j)=acos(costh)
      sint1=sin(thet1d(j))
      sint2=sin(thet2d(j))
      bu1xd(j)=sint1*cos(phi1d(j))
      bu1yd(j)=sint1*sin(phi1d(j))
      bu1zd(j)=costh
      bu2xd(j)=sint2*cos(phi2d(j))
      bu2yd(j)=sint2*sin(phi2d(j))
      bu2zd(j)=cos(thet2d(j))
      cosrtd(j)=bu1xd(j)*bu2xd(j)+bu1yd(j)*bu2yd(j)+bu1zd(j)*bu2zd(j)
      if(cosrtd(j).ge.1.0)cosrtd(j)=0.999999
      if(cosrtd(j).le.-1.0)cosrtd(j)=-0.999999
      angrtd(j)=acos(cosrtd(j))
      xsign=randproto(0)-0.5
      xsign=xsign/abs(xsign)
      if(xsign.lt.0.0)angrtd(j)=angrtd(j)-2.0*pi
      cu1xd(j)=bu1yd(j)*bu2zd(j)-bu1zd(j)*bu2yd(j)
      cu1yd(j)=-bu1xd(j)*bu2zd(j)+bu1zd(j)*bu2xd(j)
      cu1zd(j)=bu1xd(j)*bu2yd(j)-bu1yd(j)*bu2xd(j)
      ididbd(j)=1
 6230 call vecrot(bu1xd(j),bu1yd(j),bu1zd(j),cu1xd(j),cu1yd(j),cu1zd(j),
     ,  (bfracd(j)*angrtd(j)),bux,buy,buz)
      bfracd(j)=bfracd(j)+0.1
c      buang1=acos(buz)
c      buang2=atan2(buy,bux)
      if(bfracd(j).le.1.0)go to 6231
      thet1d(j)=thet2d(j)
      phi1d(j)=phi2d(j)
      bfracd(j)=0.0
      ididbd(j)=0
 6231 continue
c     Compute B field components downstream of shock in the plasma frame
c       by transforming the shock normal to the upstream plasma
c       frame, then compressing the component of B parallel to the shock
      bux=bavg*bux
      buy=bavg*buy
      buz=bavg*buz
      bup=sqrt(bux*bux+buy*buy+buz*buz)
c     Calculate upstream B field components parallel + perpendicular to shock front
c     Unit vector of shock normal at current position
      anx=cosz*cosph(j)
      any=cosz*sinph(j)
      anz=sinz
      call bdcalc(betaux(j),betauy(j),betauz(j),betau(j),gammau(j),
     ,  anx,any,anz,bux,buy,buz,bparx,bpary,bparz,bprpx,bprpy,bprpz,
     ,  bpar,bprp)
      bx(i,j)=eta*bparx+bprpx
      by(i,j)=eta*bpary+bprpy
      bz(i,j)=eta*bparz+bprpz
c      if(j.eq.12)write(6,9623)it,i,j,thet1d(j),thet2d(j),
c     ,  phi1d(j),phi2d(j),
c     ,  buang1,buang2,cosrtd(j),angrtd(j),bux,buy,buz
c     bx(i,j),by(i,j),bz(i,j)
c      if(i.eq.1.and.j.eq.12)write(6,9994)i,j,betaux(j),betauy(j),
c     ,  betauz(j),bux,bx(i,j),buy,by(i,j),buz,bz(i,j)
      n0(i,j)=eta*n0mean
c      write(6,9333)it,i,j,idelay,zcell(i,j),xcell(j),bavg,n0mean
c     Calculate the initial maximum electron energy in each cell
c     Next line relates this to direction of B field relative to shock
      gmax0(i,j)=gmaxmx*(bprp*bprp/(bpar*bpar+bprp*bprp))
      if(gmax0(i,j).lt.gmaxmn)gmax0(i,j)=gmaxmn
c     Next 3 lines assume a power-law distribution of gmax0, unrelated
c         to direction of B field
c     See http://mathworld.wolfram.com/RandomNumber.html for choosing
c      random numbers from a power-law distribution
c      xrand=randproto(0)
c      if(pexp.eq.0.0)gmax0(i,j)=gmaxmx
c      if(pexp.lt.0.0)gmax0(i,j)=((gmaxmx**pexp-gmaxmn**pexp)*xrand+
c     ,  gmaxmn**pexp)**(1.0/pexp)
      do 105 ig=1,43
  105 if(gmax0(i,j).le.gcnt(ig+1).and.gmax0(i,j).ge.gcnt(ig))
     , igcnt(ig)=igcnt(ig)+1
c
c     *** Time loop resumes here ***
c
  110 continue
c     Calculate component of magnetic field that is perpendicular to
c       the aberrated line of sight in the plasma frame
c     Line-of-sight vector in plasma frame
      slx=slos
      sly=0.0
      slz=clos
      call bcalc(betadx(i,j),betady(i,j),betadz(i,j),betad(i,j),
     ,  gammad(i,j),slx,sly,slz,bx(i,j),by(i,j),bz(i,j),dum1,dum2,
     ,  dum3,dum4,dum5,dum6,dum7,bperp(j))
      bfield(j)=sqrt(bx(i,j)**2+by(i,j)**2+bz(i,j)**2)
c      if(i.eq.1.and.j.eq.12)write(6,9996)bperp(j),bfield(j),
c     ,  bx(i,j),by(i,j),bz(i,j),
c     ,  betadx(i,j),betady(i,j),betadz(i,j),slx,sly,slz
      emisco=0.0
      ecflux=0.0
      do 111 inu=1,68
      flcomp(i,j,inu)=0.0
      flsync(i,j,inu)=0.0
      flec(i,j,inu)=0.0
      flssc(i,j,inu)=0.0
  111 flux(i,j,inu)=0.0
      bdx=betadx(i,j)
      bdy=betady(i,j)
      bdz=betadz(i,j)
      delt=(i-1)*zsize*yr*3.26/(gammad(i,j)*betad(i,j))
      bperpp=bperp(j)
      bfld=bfield(j)
c     Determine velocity of the cell plasma relative to the MD plasma
      bmparx=betamd*bdx*bdz/betad(i,j)**2
      bmpary=betamd*bdy*bdz/betad(i,j)**2
      bmparz=betamd*bdz*bdz/betad(i,j)**2
      bmprpx=-bmparx
      bmprpy=-bmpary
      bmprpz=betamd-bmparz
      betamx(j)=(bdx-bmparx-bmprpx/gammad(i,j))/(1.0d0-betamd*bdz)
      betamy(j)=(bdy-bmpary-bmprpy/gammad(i,j))/(1.0d0-betamd*bdz)
      betamz(j)=(bdz-bmparz-bmprpz/gammad(i,j))/(1.0d0-betamd*bdz)
      betamr(j)=dsqrt(betamx(j)**2+betamy(j)**2+betamz(j)**2)
      gamamr(j)=1.0d0/dsqrt(1.0d0-betamr(j)**2)
      zcl=zcell(i,j)
      zrel=zshock-0.5*zsize-zcl
      dmd(j)=sqrt(zrel**2+rcell(j)**2)
c     Calculate angle between MD seed photon and scattered l.o.s. photon
      call scatcs(bdx,bdy,bdz,clos,0.0d0,slos,xcell(j),ycell(j),
     ,   -zrel,cscat)
c     Determine Doppler factor of the MD emission in the cell's frame
      dmdp=sqrt(rcell(j)**2+(zrel/gamamr(j))**2)
      cosmd=(betamx(j)*xcell(j)+betamy(j)*ycell(j)+
     ,  betamz(j)*zrel/gamamr(j))/(betamr(j)*dmdp)
      tanmd=dsqrt(1.0d0/cosmd**2-1.0d0)
      angm=datan(tanmd)
      tanmd=dtan(2.0*datan(dtan(0.5d0*angm)/
     ,  (gamamr(j)*(1.0d0+betamr(j)))))
      cosmd=1.0d0/dsqrt(1.0d0+tanmd**2)
      deltmd(j)=1.0d0/(gamamr(j)*(1.0d0-betamr(j)*cosmd))
c     Determine time step of Mach disk seed photons from light-travel delay
      delcor=deltmd(j)/dopref
      mdmid=mdmd-(((zrf-zcl)*clos+(xrf-
     ,   xcell(j))*slos)+dmd(j))/(dtfact*zsize)
      md1=mdmid-mdrang
      if(md1.lt.1)md1=1
      if(mdmid.lt.0.or.md1.gt.mdmax)write(6,9736)it,i,j,mdmid,
     , md2,idelay,zrf,zcl,xrf,xcell(j),dmd(j),dtfact,zsize,clos,slos
      md2=mdmid+mdrang
c      write(7,9191)md1,md2,mdmid
      if(md2.gt.mdmax)md2=mdmax
      if(md1.gt.md2)md1=md2
      amdrng=md2-md1+1.0
      do 4145 inu=1,68
 4145 fmdall(inu)=0.0
      do 4147 md=md1,md2
      do 4146 inu=1,68
      syseed(inu)=0.0
      scseed(inu)=0.0
      tauxmd(inu)=0.0
 4146 ssseed(inu)=0.0
c     Mach disk's B field in frame of cell plasma transverse to cell's
c       line of sight to Mach disk (for synchrotron calculation)
      sx=rcell(j)*cosph(j)
      sy=rcell(j)*sinph(j)
      sz=zrel
      call bcalc(betamx(j),betamy(j),betamz(j),betamr(j),gamamr(j),
     ,  sx,sy,sz,bmdx(md),bmdy(md),bmdz(md),dum1,dum2,dum3,dum4,
     ,  dum5,dum6,dum7,bmperp)
c     Apply as a correction factor to previous estimate of B_perpendicular
      bpcorr=bmperp/bmdtot(md)
cc      write(5,9996)bpcorr,delcor,bmdtot(md),bmdx(md),bmdy(md),
cc     , bmdz(md),gamamr(j),betamx(j),betamy(j),betamz(j),cosbm
c
c     Calculate emission from various processes
c
c     Calculate synchrotron flux from Mach disk in frame of cell plasma
c
      nuhi=1
      do 1146 inu=1,40
      if(fsynmd(inu,md).eq.0.0.or.alphmd(inu,md).gt.9.0)
     ,   go to 1145
c     Synchrotron mean intensity for inverse Compton calculation
      ssabs=absorb(inu,md)*bpcorr**(0.5*(sen+2.0))*
     ,  delcor**(abexmd(inu,md))
      path=2.0*rsize*svmd/(tanmd*cosmd)
      if(path.gt.zsmd(md))path=zsmd(md)
      sstau=ssabs*parsec*path
      ajofnu=fsynmd(inu,md)*bpcorr**(1.0+alphmd(inu,md))*
     ,  delcor**(2.0+alphmd(inu,md))
      ajofnu=ajofnu*zsmd(md)*(emfact*parsec)
      syseed(inu)=ajofnu
      srcfn=ajofnu/sstau
      if(sstau.gt.5.0)syseed(inu)=srcfn
      if(sstau.gt.0.1.and.sstau.le.5.0)syseed(inu)=
     ,  srcfn*(1.0-exp(-sstau))
      syseed(inu)=syseed(inu)*pi*rsize**2*vmd/(dmd(j)**2)
      tauexp=1.02e4*(sen+2.0)*akapnu(nu(inu))*
     ,  bperp(j)/nu(inu)**2*parsec*dmd(j)
      tauxmd(inu)=tauexp
      if(tauexp.gt.15.0)syseed(inu)=0.0
      if(tauexp.le.15.0)syseed(inu)=syseed(inu)/exp(tauexp)
 1145 scseed(inu)=0.0
      if(inu.lt.16)fmdall(inu)=fmdall(inu)+syseed(inu)/
     ,  (amdrng/dtfact)
c      write(6,9996)nu(inu),ssabs,sstau,ajofnu,deltmd(j),n0(i,j),
c     ,  bfield(j),bpcorr,delcor,abexmd(inu,md),alphmd(inu,md),
c     ,  syseed(inu)
 1146 continue
c
c     Calculate inverse Compton flux from Mach disk in frame of cell plasma
c
      do 1147 inu=16,68
c     Inverse Compton mean intensity from Mach disk for 2nd-order
c       inverse Compton calculation in cells
      scseed(inu)=fsscmd(inu,md)*
     ,  delcor**(2.0+alfmdc(inu,md))
      scseed(inu)=scseed(inu)*zsmd(md)*pi*rsize**2*
     ,  (amjy*parsec)*vmd/dmd(j)**2
      if(tauxmd(inu).gt.15.0)scseed(inu)=0.0
      if(tauxmd(inu).le.15.0)scseed(inu)=scseed(inu)/
     ,exp(tauxmd(inu))
      fmdall(inu)=fmdall(inu)+(syseed(inu)+scseed(inu))/
     ,  (amdrng/dtfact)
 1147 continue
 4147 continue
      nuhi=0
      do 4148 inu=1,68
      ssseed(inu)=fmdall(inu)
 4148 if(fmdall(inu).gt.0.0)nuhi=inu
c     Calculate seed photon energy density energy loss calculation
      usdmd=0.0
      do 1149 inu=2,nuhi
c     Only calculate up to Klein-Nishina limit of gmin
      epslon=6.63e-27*snu(inu)/(gmin*emc2)
      if(epslon.gt.1.0)go to 1150
      if(fmdall(inu).le.SMALL_FMDALL.or.fmdall(inu-1).le.SMALL_FMDALL)go to 1148
      aaa=alog10(fmdall(inu-1)/fmdall(inu))/alog10(snu(inu)/snu(inu-1))
      usdmd=usdmd+0.5/c/(1.0-aaa)*fmdall(inu-1)*snu(inu-1)*
     ,  ((snu(inu)/snu(inu-1))**(1.0-aaa)-1.0)
 1148 continue
 1149 continue

 1150 continue
      id=(zcell(i,j)-zshock)/zsize+nzid
      if(id.lt.1)id=1
      if(id.gt.(icells+nzid))write(6,9992)i,j,id,zcell(i,j),zshock,
     ,  zsize
c     Skip flux calculation for cell if gamma_max is too low to emit at lowest frequency
c     Ratio of energy density of photons emitted by hot dust + Mach disk to
c       energy density of the magnetic field
      zdist=zdist0+(id-nend/zrat)*zsize+zshock
c     Calculate min & max angles of dust torus in plasma frame
      dphi1=asin(dtrad/sqrt(zdist**2+dtdist**2))
      dphi2=atan(dtdist/zdist)
      dth1=dphi2-dphi1
      dth2=dphi2+dphi1
      csth1=dcos(dth1)
      csth2=dcos(dth2)
      dcsth1=-(csth1-betad(i,j))/(1.0d0-csth1*betad(i,j))
      dcsth2=-(csth2-betad(i,j))/(1.0d0-csth2*betad(i,j))
      csang=0.5*(dcsth1+dcsth2)
c     Doppler factor of dust torus emission in frame of cells, used
c       to estimate frequency of peak intensity in plasma frame
      tdel=gammad(i,j)*(1.0d0-betad(i,j)*0.5d0*(csth1+csth2))
c     Correct energy density of seed photons from dust for Doppler factor
c       that include turbulent component of velocity; approximates that
c       motion is directed parallel to jet axis; also include first-order
c       dependence on scattering angle in Klein-Nishina cross-section
      useedr=useed(id)*(tdel/tdelr(id))**2/(1.0d0+betad(i,j)*csang)
      ustob=8.0*pi*(useedr+usdmd)/(bfield(j))**2
c     Calculate the maximum electron energy from gmax0 and solution to equation
c       d gamma/dt = - cc2*(b**2+8*pi*useed)*gamma**2
      tlfact=cc2*bfield(j)**2*(1.0+ustob)
      cellno=i
      tlavg=tlf1(j)/cellno
      tlf(i)=tlfact
      do 112 l=istart,i
  112 tlavg=tlavg+tlf(l)/cellno
      glim=sqrt(nu(1)*zred1/(2.8e6*bfield(j)*delta(i,j)))
      if(glim.lt.1.0)glim=1.0
      if(gammax(i-1,j).le.glim)go to 113
      tlim=(gmax0(i,j)-glim)/(tlavg*glim*gmax0(i,j))
c     skip rest of column of cells if energy losses are already severe
c      write(5,9992)it,i,j,delt,tlim,glim,gmax0(i,j),tlfact
      if(delt.lt.tlim)go to 114
  113 ididg(j)=1
c      write(5,9992)it,i,j,delt,tlim,bfield(j),bperp(j),glim,
c     ,  gmax0(i,j),gammax(i-1,j)
      go to 196
  114 if(ididg(j).eq.1)go to 196
c     iend is the last slice of cells with energetic electrons
      iend=i
c     Calculate energy distribution for the cell
      id=(zcell(i,j)-zshock)/zsize+nend
      if(id.lt.1)id=1
      if(id.gt.(icells+nend))write(6,9992)i,j,id,zcell(i,j),zshock,
     ,  zsize
      delt=zsize*yr*3.26/(gammad(i,j)*betad(i,j))
      glow=gammin(i-1,j)/(1.0+tlfact*gammin(i-1,j)*delt)
      gmrat=0.99*gammax(i-1,j)/gammin(i-1,j)
      gmratl=dlog10(gmrat)/32.0
      gmratm=dlog10(gammin(i-1,j)/glow)/11.0
      t2max=i*delt
      gamb=gammax(i-1,j)/(1.0+tlfact*gammax(i-1,j)*delt)
      ibreak=0
      if(gamb.le.glow.or.gamb.ge.gammax(i-1,j))ibreak=1
c     write(5,9994)i,j,bfield(j),n0(i,j),gmax0(i,j),usdmd,
c     ,  ustob,tlfact,delt,gamb,glow
      do 190 ie=1,44
      egam(i,j,ie)=gammin(i-1,j)*10.0**(gmratl*(ie-12))
      if(ie.lt.12)egam(i,j,ie)=glow*10.0**(gmratm*(ie-1))
      if(ie.eq.1)egam(i,j,ie)=glow+0.2*(glow*10.0**(gmratm)-glow)
      if(ibreak.eq.1.or.egam(i,j,ie).lt.gamb)go to 189
      egam(i,j,ie)=gamb
      ibreak=1
  189 ggam(ie)=egam(i,j,ie)
      enofe(i,j,ie)=0.0
      t1=(i-1)*delt
      tloss=(gmax0(i,j)-ggam(ie))/(tlavg*ggam(ie)*gmax0(i,j))
      if(tloss.le.t1)go to 188
      t2=dmin1(tloss,t2max)
      tlmin=t2max-(gmin-ggam(ie))/(tlavg*gmin*ggam(ie))
      if(ie.lt.12)t1=tlmin
      eterm1=1.0d0-ggam(ie)*t1*tlfact
      eterm2=1.0d0-ggam(ie)*t2*tlfact
      if(eterm1.lt.0.0.or.eterm2.lt.0.0)go to 188
      enofe(i,j,ie)=n0(i,j)/((sen-1.0)*tlavg*ggam(ie)**(sen+1.0))*
     ,  (eterm1**(sen-1.0)-eterm2**(sen-1.0))
c     Divide by cell crossing time since integral is over time
c      test1=eterm1**(sen-1.0)
c      test2=eterm2**(sen-1.0)
      enofe(i,j,ie)=enofe(i,j,ie)/delt
      if(enofe(i,j,ie).lt.0.0)enofe(i,j,ie)=0.0
c      if(i.eq.2)write(8,9994)j,ie,ggam(ie),enofe(i,j,ie),tloss,
c     ,  tlfact,t1,t2
  188 edist(ie)=enofe(i,j,ie)
      if((nTestOut.eq.7).and.(it.eq.1)) then
         delt_local = 0.0
         if(ie.eq.1) delt_local=delt
         write(9,14008) '188',i,j,ie,delt_local,tlfact,n0(i,j),n0mean,eta,edist(ie)
      end if
c     if(j.eq.10)write(5,9994)i,j,ggam(ie),edist(ie)
  190 continue
      gammax(i,j)=gammax(i-1,j)/(1.0d0+tlfact*delt*gammax(i-1,j))
      gammin(i,j)=gammin(i-1,j)/(1.0d0+tlfact*delt*gammin(i-1,j))
      if(gammin(i,j).lt.1.0)gammin(i,j)=1.0
      if(gammax(i,j).lt.2.0)ididg(j)=1
      if(gammax(i,j).lt.2.0)gammax(i,j)=2.0
      emold=0.0
      emeold=0.0
      emsold=0.0
c     calculate flux in mJy from cell
      ithin=0
      do 193 inu=1,22
      dustnu(inu)=dustf(id,inu)
      hnukt=4.8e-11*dustnu(inu)/tdust
      hnuktr=4.8e-11*dustnu(inu)*(tdel/tdelr(id))/tdust
      dusti(inu)=dustii(id,inu)
      if(hnuktr.gt.60.0)go to 5193
      if(dusti(inu).gt.1.0e-30)dusti(inu)=dusti(inu)*
     ,  (tdel/tdelr(id))**3*(1.0-exp(hnukt))/(1.0-exp(hnuktr))
      go to 5194
 5193 dusti(inu)=0.0
 5194 continue
  193 continue
      do 195 inu=1,68
      specin=0.0001
      restnu=nu(inu)*zred1/delta(i,j)
      snu(inu)=nu(inu)
      emisco=0.0
      fsync2(inu)=0.0
      ssabs=0.0
      ecflux=0.0
      sscflx=0.0
      if(inu.gt.44)go to 192
      emisco=ajnu(restnu)*bperp(j)*delta(i,j)**2
      fsync2(inu)=emisco*zsize*(emfact*parsec)
      fsnoab=fsync2(inu)
      if(inu.eq.1)go to 191
      if(emold.gt.0.0.and.fsync2(inu).gt.0.0)specin=
     ,  alog10(emold/fsync2(inu))/alog10(nu(inu)/nu(inu-1))
  191 continue
      if(ithin.eq.1)go to 192
      ssabs=1.02e4*(sen+2.0)*akapnu(restnu)*bperp(j)/
     ,   nu(inu)**2/delta(i,j)
      ssabs=ssabs*parsec*zsize
c     Attenuate Mach disk emission from synchrotron absorption on way to cell
      ssabsm=ssabs*dmd(j)/zsize
      if(ssabsm.ge.10.0)ssseed(inu)=0.0
      if(ssabsm.lt.10.0.and.ssabsm.gt.0.1)ssseed(inu)=
     ,     ssseed(inu)/exp(ssabsm)
c     Return to absorption within cell
c     Use rsize instead of zsize because of aberration
      ssabs=ssabs*rsize/zsize
      if(ssabs.le.(0.1/ancol))ithin=1
      srcfn=fsync2(inu)/ssabs
      if(ssabs.gt.5.0)fsync2(inu)=srcfn
      if(ssabs.gt.0.1.and.ssabs.le.5.0)fsync2(inu)=
     ,  srcfn*(1.0-exp(-ssabs))
c     Attenuation from downstream cells along l.o.s. IF significant
      tauexp=0.0
      tauexp=(nouter(j)-i)*ssabs
      if(rcell(j).gt.(0.98*rbound).and.xcell(j).le.0.0)
     ,   tauexp=0.0
      if(thlos.eq.0.0)tauexp=(nouter(j)-1)*ssabs
      if(tauexp.gt.15.0)fsync2(inu)=0.0
      if(tauexp.le.15.0)fsync2(inu)=fsync2(inu)/exp(tauexp)
      specin=alpha
  192 continue
      flsync(i,j,inu)=fsync2(inu)*(volc/zsize)*zred1/
     ,  (1.0e18*amjy*dgpc**2)*fgeom
      flux(i,j,inu)=flsync(i,j,inu)
c xxxxxxx
c      if(inu.eq.20)write(6,9994)i,j,nu(inu),flsync(i,j,inu),specin,
c     ,  bperp(j),delta(i,j),ggam(43),edist(43)
      betd=betad(i,j)
      gamd=gammad(i,j)
      if(inu.eq.1)call polcalc(bfield(j),bx(i,j),by(i,j),bz(i,j),
     ,    clos,slos,chipol)
      if(specin.lt.alpha)specin=alpha
      poldeg=(specin+1.0)/(specin+5.0/3.0)
      if(ssabs.gt.1.0)poldeg=3.0/(12.0*specin+19)
      fpol(i,j,inu)=poldeg*flsync(i,j,inu)
      pq(i,j,inu)=fpol(i,j,inu)*cos(2.0*chipol)
      pu(i,j,inu)=fpol(i,j,inu)*sin(2.0*chipol)
c      if(inu.eq.12)write(6,9992)it,i,j,bfield(j),
c     ,    bx(i,j),by(i,j),bz(i,j),clos,slos,pq(i,j,inu),
c     ,    pu(i,j,inu),fpol(i,j,inu),(chipol*rad)
      if(restnu.lt.1.0e14)go to 194
      spxec=0.0001
      spxssc=0.0001
      betd=betad(1,jcells)
      gamd=gammad(1,jcells)
      ecdust_local = ecdust(restnu)
      ecflux=ecdust_local*3.086*volc*zred1*delta(i,j)**2/dgpc**2*
     ,   fgeom
      if((nTestOut.eq.1.).and.(it.eq.1)) write(9, 9218) 192, i, j, inu, ecdust_local, delta(i,j)
c      iwrite=0
c      if(inu.eq.32)iwrite=1
c      if(inu.eq.29)write(6,9973)i,j,ecflux,delta(i,j),dusti(1),
c     ,  dusti(22),ggam(1),edist(1),ggam(44),edist(44)
      sscflx=ssc(restnu)*3.086*volc*zred1*delta(i,j)**2/dgpc**2*
     ,   fgeom
      if(nu(inu).lt.1.0e22)go to 199
      taupp=0.0
c     Pair production opacity calculation
c     Expression for anumin includes typical interaction angle
      anumin=(1.24e20/(nu(inu)*zred1))*1.24e20*(2.0*gammad(i,j))**2
      alnumn=alog10(anumin)
      inumin=(alnumn-10.0)*4+1
      if(inumin.gt.40)go to 199
      if(anumin.gt.nu(inumin))anumin=nu(inumin)
      do 197 iinu=inumin+1,40
      bp=sqrt(1.0-anumin/nu(iinu))
      if(bp.le.0.0)bp=0.001
      xsecp1=1.25e-25*(1.0-bp*bp)*((3.0-bp**4)*
     *  alog((1.0+bp)/(1.0-bp))-2.0*bp*(2.0-bp*bp))
      bp=sqrt(1.0-anumin/nu(iinu+1))
      if(bp.le.0.0)bp=0.001
      xsecp2=1.25e-25*(1.0-bp*bp)*((3.0-bp**4)*
     *  alog((1.0+bp)/(1.0-bp))-2.0*bp*(2.0-bp*bp))
      taupp=taupp+0.5*(phots(iinu)*xsecp1+phots(iinu-1)*xsecp2)*
     ,  (nu(iinu)-nu(iinu-1))*nouter(j)*zsize*parsec
  197 continue
      if(taupp.lt.0.1)go to 199
      if(taupp.lt.10.0)go to 198
      ecflux=0.0
      sscflx=0.0
      go to 199
  198 ecflux=ecflux/exp(taupp)
      sscflx=sscflx/exp(taupp)
  199 flec(i,j,inu)=ecflux
      flssc(i,j,inu)=sscflx
      flcomp(i,j,inu)=ecflux+sscflx
      flux(i,j,inu)=flsync(i,j,inu)+flcomp(i,j,inu)
      if((nTestOut.eq.1.).and.(it.eq.1)) write(9, 9218) 199, i, j, inu, ecflux, delta(i,j)
      if(emeold.gt.0.0.and.ecflux.gt.0.0)spxec=
     ,  alog10(emeold/ecflux)/alog10(nu(inu)/nu(inu-1))
      if(emsold.gt.0.0.and.sscflx.gt.0.0)spxssc=
     ,  alog10(emsold/sscflx)/alog10(nu(inu)/nu(inu-1))
  194 continue
c      if(inu.eq.21)write(5,9992)i,j,inu,ecflux,dusti(1),dusti(4),
c     ,  dusti(8),dusti(12),dusti(16),dusti(20),dusti(22)
ccc   TEST
c      if(inu.eq.25)write(5,9994)j,i,flsync(i,j,inu),bperp(j),
c     ,  n0(i,j),gammax(i,j),gammin(i,j),glow,gmin,gamb,gmax0(i,j)
      emold=fsnoab
      emeold=ecflux
      emsold=sscflx
  195 continue
cc      if(it.eq.50)write(5,9989)i,j,mdd(j),
cc     ,  n0(i,j),gammax(i,j),gammin(i,j),flsync(i,j,inu),
cc     ,  flec(i,j,52),flssc(i,j,52)
  196 continue
  200 continue
 2200 continue
c     End loop over downstream cells
c      do 7999 ig=1,44
c 7999 write(5,9993)gcnt(ig),igcnt(ig)
  299 tflold=0.0
      do 500 inu=1,68
      tflux=0.0
      tsflux=0.0
      tcflux=0.0
      tecfl=0.0
      tsscfl=0.0
      alph=0.0
      qcum=0.0
      ucum=0.0
      iwp=0
      if(it.eq.1)iwp=1
      if(it.eq.25.or.it.eq.50)iwp=1
      if(it.eq.75.or.it.eq.100)iwp=1
      if(it.eq.125.or.it.eq.150)iwp=1
      if(it.eq.175.or.it.eq.200)iwp=1
      if(it.eq.225.or.it.eq.250)iwp=1
      if(it.eq.275.or.it.eq.300)iwp=1
      if(it.eq.325.or.it.eq.350)iwp=1
      if(it.eq.375.or.it.eq.400)iwp=1
      if(it.eq.425.or.it.eq.450)iwp=1
      if(it.eq.475.or.it.eq.500)iwp=1
      if(it.eq.525.or.it.eq.550)iwp=1
      if(it.eq.575.or.it.eq.600)iwp=1
      if(it.eq.625.or.it.eq.650)iwp=1
      if(it.eq.675.or.it.eq.700)iwp=1
      if(it.eq.725.or.it.eq.750)iwp=1
      if(it.eq.775.or.it.eq.800)iwp=1
      if(it.eq.825.or.it.eq.850)iwp=1
      if(it.eq.875.or.it.eq.900)iwp=1
      if(it.eq.925.or.it.eq.950)iwp=1
      if(it.eq.975.or.it.eq.1000)iwp=1
      if(iwp.ne.1.or.inu.ne.3)go to 7300
      write(filnam,"(a7,i4.4,a4)") "temzmap",it,".txt"
      open (7,iostat=ios, err=9000, file=dpath//filnam,status=fstat)
      write(7,9875)
 7300 continue
      do 300 j=1,(jcells-1)
      do 300 i=1,icelmx(j)
c      i=1
      qcum=qcum+pq(i,j,inu)
      ucum=ucum+pu(i,j,inu)
      tflux=tflux+flux(i,j,inu)
      tsflux=tsflux+flsync(i,j,inu)
      if(nTestOut.eq.5) write(9,14005) it,i,j,0,inu,'flsync',flsync(i,j,inu),'tsflux',tsflux
      tcflux=tcflux+flcomp(i,j,inu)
      tecfl=tecfl+flec(i,j,inu)
      tsscfl=tsscfl+flssc(i,j,inu)
      if(inu.ne.3)go to 298
c     Position of cell on sky, in milliarcseconds (for z=0.859)
      xobs=(zcell(i,j)*slos+xcell(j)*clos)/7.7
      yobs=ycell(j)/7.7
      poldgg=0.0
      pang=0.0
      if(flsync(i,j,inu).le.1.0e-6)go to 298
      if(iwp.ne.1)go to 298
      poldg=100.0*fpol(i,j,inu)/flsync(i,j,inu)
      pangg=0.5*atan2(pu(i,j,inu),pq(i,j,inu))*rad
      write(7,9876)i,j,xobs,yobs,(0.001*flsync(i,j,inu)),
     ,  (0.001*pq(i,j,inu)),(0.001*pu(i,j,inu)),poldg,pangg,
     ,  tflux,qcum,ucum,betad(i,j),it
c      if(inu.eq.9)write(5,9994)j,i,flux(i,j,inu),flsync(i,j,inu),
c     ,  flcomp(i,j,inu),flec(i,j,inu),flssc(i,j,inu)
  298 continue
  300 continue
c      phots(inu)=tflux*amjy*
c     ,  (4*pi*dgpc**2*(zsize*1.0e18))/
c     , ((6.63e-27*nu(inu))*c*volc*jcells*zred1)
      phots(inu)=0.0
      phalph(inu)=alph+1.0
      if(tflold.eq.0.0.or.tflux.eq.0.0)go to 390
      if(inu.gt.1)alph=alog10(tflold/tflux)/
     ,alog10(nu(inu)/nu(inu-1))
      poldeg=0.0
      if(tflux.gt.0.0)poldeg=sqrt(qcum**2+ucum**2)/tflux
      pang=0.5*atan2(ucum,qcum)*rad
      tflux=0.001*nu(inu)*tflux
      tsflux=0.001*nu(inu)*tsflux
      tcflux=0.001*nu(inu)*tcflux
      tecfl=0.001*nu(inu)*tecfl
      tsscfl=0.001*nu(inu)*tsscfl
      tsscf2=0.001*nu(inu)*tsscf2
c     Print out SED for every ispecs-th time step
  390 continue
c  390 if(ispec.gt.1)go to 391
      timeo=it*dtime
      if(inu.eq.1)write(3,9990)timeo
c     Write SED to file
      tfsync=1.0e3*tsflux/nu(inu)
      if(inu.eq.20)tfsyno=1.0e3*tsflux/nu(inu)
      write(3,9996)nu(inu),tflux,alph,tfsync,tecfl,tsscfl
  391 if(inu.eq.19)tfl19=tflux
      if(inu.eq.20)tfl20=tflux
      if(inu.eq.21)tfl21=tflux
      if(inu.eq.32)tfl32=tflux
      if(inu.eq.33)tfl33=tflux
      if(inu.eq.34)tfl34=tflux
      if(inu.eq.33)tfcomx=1.0e6*tcflux/nu(inu)
      if(inu.eq.33)tfec=tecfl
      if(inu.eq.33)tfssc=tsscfl
      if(inu.eq.52)tfl52=tflux
      if(inu.eq.53)tfl53=tflux
      if(inu.eq.54)tfl54=tflux
      if(inu.eq.53)tfcomp=1.0e9*tcflux/nu(inu)
      if(inu.eq.53)tfec=tecfl
      if(inu.eq.53)tfssc=tsscfl
      if(inu.eq.11)tfl11=tflux
      if(inu.eq.6)tfl6=tflux
      tflold=1000.0*tflux/nu(inu)
      if(inu.eq.3)pdeg3=100.0*poldeg
      if(inu.eq.3)pang3=pang
      if(inu.eq.6)pdeg6=100.0*poldeg
      if(inu.eq.6)pang6=pang
      if(inu.eq.8)pdeg8=100.0*poldeg
      if(inu.eq.8)pang8=pang
      if(inu.eq.12)pdeg12=100.0*poldeg
      if(inu.eq.12)pang12=pang
      if(inu.eq.16)pdeg16=100.0*poldeg
      if(inu.eq.16)pang16=pang
      if(inu.eq.20)pdeg20=100.0*poldeg
      if(inu.eq.20)pang20=pang
      if(inu.eq.24)pdeg24=100.0*poldeg
      if(inu.eq.24)pang24=pang
      if(inu.eq.28)pdeg28=100.0*poldeg
      if(inu.eq.28)pang28=pang
      if(inu.eq.32)pdeg32=100.0*poldeg
      if(inu.eq.32)pang32=pang
  500 continue
      close(7)
      alp20a=0.0
      alp20b=0.0
      alp33a=0.0
      alp33b=0.0
      alp53a=0.0
      alp53b=0.0
      if(tfl19.gt.0.0.and.tfl20.gt.0.0)alp20a=alog10(tfl19/tfl20)/
     ,alog10(nu(20)/nu(19))+1.0
      if(tfl20.gt.0.0.and.tfl21.gt.0.0)alp20b=alog10(tfl20/tfl21)/
     ,alog10(nu(21)/nu(20))+1.0
      if(tfl32.gt.0.0.and.tfl33.gt.0.0)alp33a=alog10(tfl32/tfl33)/
     ,alog10(nu(33)/nu(32))+1.0
      if(tfl33.gt.0.0.and.tfl34.gt.0.0)alp33b=alog10(tfl33/tfl34)/
     ,alog10(nu(34)/nu(33))+1.0
      if(tfl52.gt.0.0.and.tfl53.gt.0.0)alp53a=alog10(tfl52/tfl53)/
     ,alog10(nu(53)/nu(52))+1.0
      if(tfl53.gt.0.0.and.tfl54.gt.0.0)alp53b=alog10(tfl53/tfl54)/
     ,alog10(nu(54)/nu(53))+1.0
      alph20=0.5*(alp20a+alp20b)
      alph33=0.5*(alp33a+alp33b)
      alph53=0.5*(alp53a+alp53b)
      timeo=it*dtime
c     Write light curve points to file
      write(4,9991)it,timeo,nu(20),tfl20,alph20,tfsyno,nu(33),
     ,  tfl33,alph33,tfcomx,nu(53),tfl53,alph53,tfcomp,tfec,
     ,  tfssc,nu(6),tfl6,ncells
c     Write selected polarization data to file
      write(5,9988)it,timeo,nu(3),pdeg3,pang3,nu(6),pdeg6,pang6,
     ,  nu(8),pdeg8,pang8,nu(12),pdeg12,pang12,
     ,  nu(16),pdeg16,pang16,nu(20),pdeg20,pang20,
     ,  nu(24),pdeg24,pang24
      print *, 'Done with time loop', it
c     Set up next time step by shifting physical conditions 1 slice down jet
      if(it.eq.itlast)go to 9000
      do 598 j=1,jcells
      do 598 i=icelmx(j)+1,istart,-1
      gmax0(i,j)=gmax0(i-1,j)
      bx(i,j)=bx(i-1,j)
      by(i,j)=by(i-1,j)
      bz(i,j)=bz(i-1,j)
      n0(i,j)=n0(i-1,j)
      betadx(i,j)=betadx(i-1,j)
      betady(i,j)=betady(i-1,j)
      betadz(i,j)=betadz(i-1,j)
      betad(i,j)=betad(i-1,j)
      gammad(i,j)=gammad(i-1,j)
      delta(i,j)=delta(i-1,j)
  598 continue
c      write(3,8889)
      ispec=ispec-1
      if(ispec.eq.0)ispec=ispecs+1
c     Move cells in time array to make room for next time step
      do 599 inu=1,68
      do 599 md=1,(mdmax-1)
      alphmd(inu,md)=alphmd(inu,md+1)
      alfmdc(inu,md)=alfmdc(inu,md+1)
      fsynmd(inu,md)=fsynmd(inu,md+1)
      fsscmd(inu,md)=fsscmd(inu,md+1)
      absorb(inu,md)=absorb(inu,md+1)
  599 continue
      go to 9
 6666 format('# freq(Hz)  flux density (Jy Hz) spectral index')
 6667 format('#'/'#   time(d)     freq(Hz)  F(Jy Hz)  alpha ',
     ,   '    F(mJy)   ',
     ,   ' freq(Hz)  F(Jy Hz)  alpha     F(muJy) ',
     ,   '  freq(Hz)  F(Jy Hz)   alpha     F(nJy)    F(EC)   ',
     ,      '  F(SSC)  no. of live cells')
 6668 format('#'/'#     time(d)  freq(Hz)   p(%)   chi(deg) ',
     ,   '  freq      p      chi     freq       p      chi  ',
     ,   '   freq      p       chi     freq       p      chi  ',
     ,   '    freq       p      chi    freq       p       chi')
 7777 format(i4,1x,i4,2f12.4,2f8.2,2x,f10.4)
 8888 format(i5,f9.2)
 8889 format(/)
 9111 format(a10/i3/f5.3/f5.3/f4.2/f4.2/f5.3/f3.1/e7.1/f5.3/f7.1/
     ,f5.1/f6.1/d9.5/d9.5/d6.1/d6.1/d6.2/f6.1/f4.1/f3.1/f3.1/f3.1/e7.2)
 9215 format('md',i6,' inu',i5, ' fsscmd(inu,md)',es12.5)
 9216 format('md',i6,' inu',i5, ' fsynmd(inu,md)',f15.2)
 9217 format('md',i6,' inu',i5, ' fsync(inu)',f15.3)
 9218 format(i4,' i ', i5, ' j ', i5, ' inu ', i5, ' ecflux ', es12.4,' delta ', es12.4)
 9219 format('j ', i5, ' imax ', i5)
 9220 format('ncells ', i8, ' j ', i5)
 9221 format('1908 ncells ', i8, ' j ', i5, ' i ', i5)
 9222 format('idelay out of bounds ',2i12,3i6,1p7e11.3)
 9223 format('*idelay out of bounds ',2i12,3i6,1p7e11.3)
 9232 format(2i4,3i6,1p12e10.2)
 9299 format(i12,2i5,2i6,3i10,1p5e12.4)
 9333 format(4i7,1p12e10.2)
 9623 format('* ',3i5,13f8.4)
 9736 format(6i7,1p10e12.3)
 9875 format('#    i     j    x(mas)     y(mas)   flux(Jy)     Q(Jy)',
     ,   '       U(Jy)      P(%)    EVPA(deg)',
     ,   '  tot flux    qcum      ucum     betad      it')
 9876 format(2i5,1x,1p11e11.3,i5)
 9891 format('Program halted: betadd < sound speed ',1p6e12.5)
 9911 format('seed:  ',2i5,1p11e10.3)
 9988 format(i5,f8.2,2x,7(1pe8.2,1x,0pf7.3,1x,f8.3,1x))
 9989 format(i6,2x,i6,2x,i6,1p16e10.2)
 9990 format(///'Time = ',f8.2,' days'/'  freq(Hz))',2x,
     , 'Ftot(Jy Hz)  sp. index    Fsynch',4x,'     F(EC)',
     , 4x,'F(SSC-MD)')
 9991 format(i5,f8.2,2x,1p16e10.2,1x,i7)
 9992 format('** ',3i6,1p12e12.4)
 9993 format(e12.4,1x,i10,1x,2f8.3)
 9994 format(i5,2x,i5,1p16e12.4)
 9995 format('E dist',2x,1p2e11.4)
 9996 format(1p12e12.4)
 9997 format(3e10.2,2x,1p10e12.4)
 9998 format('# mean, std. dev. pol and EVPA',4f12.5/
     , '# mean, std. dev. pol., layers 1-7',2f12.5)
 9999 format(10f10.3)
14001 format(i5,' i',i5,' j',i5,' md',i6,' inu',i5,' ',a,f10.4)
14005 format(i5,' i',i5,' j',i5,' md',i6,' inu',i5,' ',a,f10.4,' ',a,f10.4)
14006 format(i7,' ',i8)
14007 format(i7,' ',es12.4)
14008 format(a5,' ',i7,' ',i7,' ',i2,' ',i8,' ',es12.4,' ',es12.4,' ',es12.3,' ',es12.3,' ',es12.4,' ',es12.5)
14009 format(a10,' ',es12.4)
      close (4, status='keep')
      close (3, status='keep')
      close (9, status='keep')
 9000 stop
      end
c
c     Subroutine that calculates the seed photon intensity from emission by hot dust
      function seedph(f)
      common/cdust/csth1,csth2,dcsth1,dcsth2,dsang,tdust
      common/cvel/bdx,bdy,bdz,gammad,betad
      common/cfreq/freq
      real*8 bdx,bdy,bdz,gammad,betad,csth1,csth2,
     ,  dcsth1,dcsth2,cs1,cs2,csm
      external sdgran
      freq=f
      cs1=csth2
      cs2=csth1
      csm=0.5d0*(cs1+cs2)
      call qg5(cs1,csm,sdgran,ans1)
      call qg5(csm,cs2,sdgran,ans2)
c     Multiply by 2pi, but then divide by 4pi to get mean intensity
      seedph=ans1+ans2
      return
      end
c
      function sdgran(cs)
      common/cdust/csth1,csth2,dcsth1,dcsth2,dsang,tdust
      common/cvel/bdx,bdy,bdz,gammad,betad
      common/cfreq/freq
      real*8 cs,csp,bdx,bdy,bdz,gammad,betad,csth1,csth2,
     ,  dcsth1,dcsth2
      hok=4.80e-11
      sdgran=0.0
      csp=-(cs-betad)/(1.0d0-betad*cs)
      tdel=1.0/(gammad*(1.0d0-betad*csp))
      f=freq/tdel
      expon=hok*f/tdust
c      write(6,9000)cs,tdel,freq,f,expon
      if(expon.lt.0.01)go to 10
      if(expon.gt.10.0)go to 20
      val=(1.33e-26*f)*f/((9.0e20/f)*(exp(expon)-1.0))
      go to 25
   10 val=(2.76e-16*f)*tdust*(f/9.0e20)
      go to 25
   20 if(expon.gt.20.0)go to 26
      val=(1.33e-26*f)*f/((9.0e20/f)*exp(expon))
   25 sdgran=val*tdel*tdel
   26 return
c 9000 format('++ ',1p10e12.4)
      end
c
c     Subroutine to calculate synchrotron emission coefficient
      function ajnu(anu)
      common/cparm/zred1,bfield,b
      common/cdist/gam,edist
      dimension gam(44),edist(44)
      real*8 gam
      ajnu=0.0
      if(b.lt.1.0e-5)go to 1000
      xfact=2.38e-7*anu/b
      x=xfact/(gam(1)*gam(1))
      if(x.gt.0.01)go to 10
      rfact=1.5*x**0.33333333
      go to 20
   10 rfact=0.0
      if(x.lt.25.0)rfact=(1.08895*x**0.20949-
     , 2.35861e-3/x**0.79051)/exp(x)
   20 gran1=edist(1)*rfact
      do 100 i=2,44
      x=xfact/(gam(i)*gam(i))
      if(x.gt.0.01)go to 30
      rfact=1.5*x**0.33333333
      go to 40
   30 rfact=0.0
      if(x.lt.25.0)rfact=(1.08895*x**0.20949-
     , 2.35861e-3/x**0.79051)/exp(x)
   40 gran2=edist(i)*rfact
      if(gran1.eq.0.0.or.gran2.eq.0.0)go to 45
      grat=gam(i)/gam(i-1)
      alg=alog10(grat)
      if(alg.eq.0.0)go to 45
      a=1.0+alog10(gran2/gran1)/alg
      if(a.gt.5.0.or.a.lt.-5.0)go to 45
      if(a.lt.0.01)go to 45
      addit=gran1*(grat**a-1.0)*gam(i-1)/a
      go to 46
   45 addit=0.5*(gran1+gran2)*(gam(i)-gam(i-1))
   46 ajnu=ajnu+addit
c      if(anu.gt.9.0e11.and.anu.lt.1.1e12)write(5,1111)
c     , anu,b,x,gam(i),edist(i),gran1,gran2,a,addit,ajnu,rfact
      gran1=gran2
  100 continue
 1111 format(1p14e11.4)
 1000 return
      end
c
c     Subroutine to calculate synchrotron absorpion coefficient
      function akapnu(anu)
      common/cparm/zred1,bfield,b
      common/cdist/gam,edist
      dimension gam(44),edist(44)
      real*8 gam
      akapnu=0.0
      xfact=2.38e-7*anu/b
      x=xfact/(gam(1)*gam(1))
      akapnu=0.0
      f1=0.0
      if(x.gt.25.0)go to 10
      f1=1.8*x**0.3/exp(x)
   10 gran1=f1*edist(1)/gam(1)
      do 100 i=2,44
      x=xfact/(gam(i)*gam(i))
      gran2=0.0
      if(x.gt.25.0)go to 45
      f2=1.8*x**0.3/exp(x)
      gran2=f2*edist(i)/gam(i)
      if(gran1.eq.0.0.or.gran2.eq.0.0)go to 45
      grat=gam(i)/gam(i-1)
      alg=alog10(grat)
      if(alg.eq.0.0)go to 45
      a=1.0+alog10(gran2/gran1)/alg
      if(a.gt.5.0.or.a.lt.-5.0)go to 45
      if(a.lt.0.01)go to 45
      addit=gran1*(grat**a-1.0)*gam(i-1)/a
      go to 46
   45 addit=0.5*(gran1+gran2)*(gam(i)-gam(i-1))
   46 akapnu=akapnu+addit
c      if(anu.gt.9.1e10.and.anu.lt.1.1e11)write(3,1111)
c     , anu,b,x,gam(i),edist(i),gran1,gran2,a,addit,ajnu,rfact
      gran1=gran2
  100 continue
 1111 format(1p14e11.4)
 1000 return
      end
c     Subroutine to carry out 10-point Gaussian integration
      subroutine qg10(a,b,func,ss)
      real ss,func
      real*8 a,b
      external func
c     Returns as ss the integral of the function func between a and 
c     b, by ten-point Gauss-Legendre integration:
c     the function is evaluated exactly ten
c     times at interior points in the range of integration.
      integer j
c     The abscissas and weights:
      REAL dx,xm,xr,w(5),x(5) 
      SAVE w,x
      DATA w/.2955242247,.2692667193,.2190863625,.1494513491,
     ,  .0666713443/
      DATA x/.1488743389,.4333953941,.6794095682,.8650633666,
     ,  .9739065285/
      xm=0.5*(b+a)
      xr=0.5*(b-a)
      ss=0
c     Will be twice the average value of the function, since the ten
c     weights (five numbers above each used twice) sum to 2. 
      do 11 j=1,5
      dx=xr*x(j)
      ss=ss+w(j)*(func(xm+dx)+func(xm-dx))
   11 continue
c     Scale the answer to the range of integration.
      ss=xr*ss 
      return
      END
c     Subroutine to carry out 5-point Gaussian integration
      subroutine qg5(a,b,func,ss)
      real*4 ss,func
      external func
c     Returns as ss the integral of the function func between a and 
c     b, by 5-point Gauss-Legendre integration:
c     the function is evaluated exactly ten
c     times at interior points in the range of integration.
      integer j
c     The abscissas and weights:
      real*8 a,b,dx,xm,xr,w(5),x(5) 
      save w,x
      data w/.236926885,.478628670,.568888889,.47862867,
     ,  .236926885/
      data x/-.906179846,-.538469310,0.0,.53846931,.906179846/
      ss=0
      xm=0.5*(a+b)
      xr=0.5*(b-a)
      do 11 j=1,5
      ss=ss+w(j)*func(xm+x(j)*xr)
   11 continue
c     Scale the answer to the range of integration.
      ss=xr*ss 
      return
      end
c     Subroutine to carry out 3-point Gaussian integration     
      subroutine qg3(xl,xu,fct,y)
c     Gaussian quadrature integration using 3 points
      external fct
      a=0.5*(xu+xl)
      b=xu-xl
      y=0.3872983*b
      y=0.2777778*(fct(a+y)+fct(a-y))
      y=b*(y+0.4444444*fct(a))
      return
      end
c     The next 4 subroutines are from Ritaban Chatterjee to create variations
c       according to an input PSD with slopes beta1 and beta2 at
c       variational frequencies below and above a break frequency nu_break
C**************************************************************
      subroutine psdsim(N,beta1,beta2,nu_break,t_incre1,lc_sim)
C**************************************************************
C N=number of data points in the lc, should be an integer power of 2, N<=32768
C t_incre1=increment in time in each step while resampling the simulated data. 
C This must be larger than the smallest interval between successive data points in the input light curve

      implicit none

      real*8 nu(131072),flux_s1(131072),dat(262144),R(262144)
      real*8 dataim_s1(131072),datareal_s1(131072),amp(131072)
      real*8 flux_s2(131072),dataim_s2(131072),datareal_s2(131072)
      real*8 dataim(131072),datareal(131072),flux_s(131072)
      real*8 fac_norm,fac_norm2
      real*4 ann,beta1,beta2,nu_break,t_incre1,lc_sim(131072)
      integer N,j,i,ifile,ik,nn,ISEED1,ISEED2,isign

      fac_norm=1./(N*t_incre1)
      fac_norm2=N**2./(2.*N*t_incre1)

      ISEED1=393521
      ISEED2=17263
      call RNE2IN(ISEED1,ISEED2)

        call RNSTNR(R,N/2)
        call RNE2OT(ISEED1,ISEED2)
        do j=1,N/2  
            nu(j)= j*fac_norm/86400.                      
            if (nu(j) .le. nu_break) then
              flux_s(j)=nu(j)**beta1
              datareal(j)=dsqrt(fac_norm2*0.5*flux_s(j))*R(j)
            else
              flux_s(j)=(nu(j)**beta2)
     >                      *((nu_break**beta1)/(nu_break**beta2))
              datareal(j)=dsqrt(fac_norm2*0.5*flux_s(j))*R(j)
            endif
        end do
        call RNSTNR(R,N/2)
        call RNE2OT(ISEED1,ISEED2)

        do j=1,N/2
            dataim(j)=dsqrt(0.5*fac_norm2*flux_s(j))*R(j)
        end do

           dat(1)=0.0
           dat(2)=0.0
           dat(N+1)=datareal(N/2)
           dat(N+2)=0 ! this is at +/-1/T and has to obey A+iB = A-iB => B=0          
           j=1

        do i=3,N-1,2
             dat(i)=datareal(j)
             dat(2*(N+1)-i)=dat(i)
             dat(i+1)=dataim(j)
             dat(2*(N+1)-i+1)=-1.0*dat(i+1)
           j=j+1
        end do
        nn=N
        isign=-1
        call four1(dat,nn,isign)
        j=1
        do i=1,2*nn,2
           ann=nn
           lc_sim(j)=dat(i)/ann !taking the real part
C           write(8,*)lc_sim(j)
           j=j+1
        end do

      return
      end  
C*******************************************
      SUBROUTINE RNSTNR(R,N)
C*******************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(N), U(2)
      I=1
10    CONTINUE
      CALL RNECUY(U,2)
      V1=2.*U(1)-1.
      V2=2.*U(2)-1.
      S=V1*V1+V2*V2
      IF(S.GE.1.) GO TO 10
      ROOT=SQRT(-2.*LOG(S)/S)
      R(I)=V1*ROOT
      IF(I.LT.N) R(I+1)=V2*ROOT
      I=I+2
      IF(I.LE.N) GO TO 10
      END

C*******************************************     
      SUBROUTINE RNECUY(U,N)
C*******************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(N)
      PARAMETER (M1=2147483563,M1MIN1=2147483562,IA1=40014,
     +           IQ1=53668,IR1=12211)
      PARAMETER (M2=2147483399,IA2=40692,IQ2=52774,IR2=2791)
      PARAMETER (XM1INV=4.656613D-10)
      DATA IX1,IX2 /123456, 654321/
      DO 10 I=1,N
C Produce integer random number X1 from first MLCG
        K=IX1/IQ1
        IX1=IA1*(IX1-K*IQ1)-K*IR1
        IF (IX1.LT.0) IX1=IX1+M1
C Produce integer random number X2 from second MLCG
        K=IX2/IQ2
        IX2=IA2*(IX2-K*IQ2)-K*IR2
        IF (IX2.LT.0) IX2=IX2+M2
C Combine
        IZ=IX1-IX2
        IF(IZ.LT.1) IZ=IZ+M1MIN1
C Normalize and transform to floating point
        U(I)=DBLE(IZ)*XM1INV
   10 CONTINUE
      RETURN
C Input of seed
      ENTRY RNE2IN(ISEED1,ISEED2)
      IX1=ISEED1
      IX2=ISEED2
      RETURN
C Output of seed
      ENTRY RNE2OT(ISEED1,ISEED2)
      ISEED1=IX1
      ISEED2=IX2
      RETURN
      END
                                                           
C*****************************************************************
       SUBROUTINE four1(data,nn,isign)
C****************************************************************

       INTEGER isign,nn
       REAL*8 data(2*nn)
c  Replaces data(1:2*nn) by its discrete Fourier transform, 
c  if isign is input as 1; or replaces
c  data(1:2*nn) by nn times its inverse discrete Fourier 
c  transform, if isign is input as -1.
c  data is a complex array of length nn or, equivalently, 
c  a real array of length 2*nn.
c  nn MUST be an integer power of 2 (this is not checked!).
       INTEGER i,istep,j,m,mmax,n
       REAL*8 tempi,tempr
       real*8 theta,wi,wpi,wpr,wr,wtemp !Double precision 
                              !   for the trigonometric recurrences. 
       n=2*nn
       j=1
       do i=1,n,2 !This is the bit-reversal section of the routine.
         if(j.gt.i)then
            tempr=data(j) ! Exchange the two complex numbers.
            tempi=data(j+1)
            data(j)=data(i)
            data(j+1)=data(i+1)
            data(i)=tempr
            data(i+1)=tempi
         endif
       m=nn
1      if ((m.ge.2).and.(j.gt.m)) then
         j=j-m
         m=m/2
         goto 1
       endif
         j=j+m
       end do 
       mmax=2 ! Here begins the Danielson-Lanczos section of 
              ! the routine.
2      if (n.gt.mmax) then ! Outer loop executed log2 nn times.
           istep=2*mmax
           theta=6.28318530717959d0/(isign*mmax) !Initialize for 
                                  ! the trigonometric recurrence.
           wpr=-2.d0*sin(0.5d0*theta)**2
           wpi=sin(theta)
           wr=1.d0
           wi=0.d0
           do m=1,mmax,2 ! Here are the two nested inner loops. 

              do i=m,n,istep
                j=i+mmax        !This is the Danielson-Lanczos formula:
                tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
                tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
                data(j)=data(i)-tempr
                data(j+1)=data(i+1)-tempi
                data(i)=data(i)+tempr
                data(i+1)=data(i+1)+tempi
              end do
            wtemp=wr                    !Trigonometric recurrence.
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
          end do
        mmax=istep
        goto 2                          !Not yet done.
      endif                             !All done.
      return
      END
c
c     polcalc computes the polarization angle chi based on relations
c       derived by Lyutikov et al. 2005, MNRAS, 360, 869
c
      Subroutine polcalc(b,bx,by,bz,clos,slos,chi)
      common/cvel/bdx,bdy,bdz,gammad,betad
      real*8 clos,slos,bdx,bdy,bdz,gammad,betad,term,ndq,q2,
     ,  qx,qy,qz,ex,ey,ez,term1,term2,term3
      bxh=bx/b
      byh=by/b
      bzh=bz/b
      term1=slos*bxh+clos*bzh
      term2=slos*bdx+clos*bdz
      term3=(bdx*bxh+bdy*byh+bdz*bzh)*gammad/(1.0d0+gammad)
      qx=bxh+(term1-term3)*bdx-term2*bxh
      qy=byh+(term1-term3)*bdy-term2*byh
      qz=bzh+(term1-term3)*bdz-term2*bzh
      q2=qx*qx+qy*qy+qz*qz
      ndq=qx*slos+qz*clos
      term1=-qy*clos
      term2=qx*clos-qz*slos
      term3=qy*slos
      term=dsqrt(q2-ndq*ndq)
c      ex=term1/term
      ey=term2/term
c      ez=term3/term
c      term=ex*ex+ey*ey+ez*ez
c      coschi=(-clos*ex+slos*ez)
c      chi=acos(coschi)
c      term2=ey/dsqrt(term)
      if(ey.gt.1.0d0.and.ey.lt.1.000001d0)ey=0.99999999d0
      if(ey.lt.-1.0d0.and.ey.gt.-1.000001d0)ey=-0.99999999d0
      chi=asin(ey)
      return
      end
      ! Function to call rand() most of the time, but when fixedRandFlag is non-zero, 
      ! get the random values from a file. When all the numbers are used up, start
      ! back at the beginning of the number list
      function randproto(seed)
      common/cfixedrand/fixedRandFlag, fixedRandFileOpened, fixedRandData, fixedRandCounter
      dimension fixedRandData(100000)
      real fixedRandData
      integer fixedRandFlag, fixedRandFileOpened, fixedRandCounter
      integer*4 seed

      if(fixedRandFileOpened.eq.0) then
         ! need to open the file with all the pre-generated random nums
         ! and load into a shared array
         ! At some point, I want to have the ability to set an environment variable for the filename:
         !character(len=512) :: filename
         !call get_environment_variable('BLZFIXEDRAND', filename, 512, status)
         !if(status.eq.1) then
         !   filename = './fixedrand.dat'
         !endif
         open(99, FILE = 'fixedrand.dat') ! At the moment, as you can see, this file must be in running/current directory
         do i=1,100000
            read(99, *) fixedRandData(i)
         end do
         fixedRandFileOpened = 1
         fixedRandCounter = 1
         close(99)
      endif

      if(fixedRandFlag.eq.0) then
         randproto = rand(seed)
      else
         randproto = fixedRandData(fixedRandCounter)
         if(fixedRandCounter.ge.100000) then
            fixedRandCounter = 1
         else
            fixedRandCounter = fixedRandCounter + 1
         endif
      endif
      return
      end
c
c     vdcalc computes downstream velocity vector; Lorentz transformation of unit
c     vector along shock front follows Lyutikov et al. (2003, ApJ, 597, 998)
c
      subroutine vdcalc(vx,vy,vz,sx,sy,sz,vdx,vdy,vdz,vd,gd,eta)
      common/cang/cosz,sinz
      real*8 vx,vy,vz,sx,sy,sz,v2,v,g,spx,spy,spz,s,g2,
     ,  dotprd,vparx,vpary,vparz,vprpx,vprpy,vprpz,vprp2,vprp,uprp,
     ,  vd,gd,vd2,vdx,vdy,vdz,vfact,vdprpx,vdprpy,vdprpz,vdprp2,
     ,  gdprp,eta,cosz,sinz
      v2=vx*vx+vy*vy+vz*vz
      v=dsqrt(v2)
      g2=1.0d0/(1.0d0-v2)
      g=dsqrt(g2)
      dotprd=vx*sx+vy*sy+vz*sz
c     Calculate components of velocity vector parallel and perpendicular to shock front
      vparx=dotprd*sx
      vpary=dotprd*sy
      vparz=dotprd*sz
      vprpx=vx-vparx
      vprpy=vy-vpary
      vprpz=vz-vparz
      vprp2=vprpx*vprpx+vprpy*vprpy+vprpz*vprpz
      vprp=dsqrt(vprp2)
c     Shock jump condition, from Konigl (1980, Phys. Fluids, 23, 1083)
      vfact=1.0d0
      uprp=g*vprp
c     uprp must exceed the proper sound speed, 1/sqrt(2) for a shock
c     Otherwise, it is a sound wave and the velocity does not change significantly
      if(uprp.gt.0.7071e0)vfact=(1.0d0+1.0d0/(g2*vprp2))/(3.0d0)
      vdprpx=vprpx*vfact
      vdprpy=vprpy*vfact
      vdprpz=vprpz*vfact
      vdprp2=vdprpx*vdprpx+vdprpy*vdprpy+vdprpz*vdprpz
c     Shock compression ratio (downstream density/upstream density)
c       From Cawthorne + Cobb (1990, ApJ, 350, 536)
      eta=dsqrt((8.0d0*v2*sinz*sinz-1.0d0/g2)/(1.0d0-v2*cosz*cosz))*
     ,  v*g*sinz
c     eta must be >= 1; if formula gives < 1, it is not a shock
      if(eta.lt.1.0)eta=1.0
c      eta=dsqrt(vprp2*(1.0d0-vdprp2)/(vdprp2*(1.0d0-vprp2)))
      vdx=vparx+vdprpx
      vdy=vpary+vdprpy
      vdz=vparz+vdprpz
      vd2=vdx*vdx+vdy*vdy+vdz*vdz
      vd=dsqrt(vd2)
      gd=1.0d0/dsqrt(1.0d0-vd2)
c      write(5,9999)v,g,vd,gd,eta
c 9999 format(1p5e12.4)
      return
      end
c
c     bdcalc computes magnetic field components downstream of shock
c
      subroutine bdcalc(vx,vy,vz,v,g,ax,ay,az,bx,by,bz,
     ,  bparx,bpary,bparz,bprpx,bprpy,bprpz,bpar,bprp)
      real*8 vx,vy,vz,ax,ay,az,sx,sy,sz,a,v2,v,g,g1,g2,gp1,
     ,  apx,apy,apz
c     v is velocity in units of c, a is a unit vector normal to the shock front
c     b is magnetic field vector
      v2=vx*vx+vy*vy+vz*vz
      v=dsqrt(v2)
      g2=1.0d0/(1.0d0-v2)
      g=dsqrt(g2)
      g1=g-1.0d0
      gp1=g+1.0d0
c     Determine unit vector of shock normal. in plasma frame
      dotprd=vx*ax+vy*ay+vz*az
      apx=g*ax-g1*dotprd*vx/v2
      apy=g*ay-g1*dotprd*vy/v2
      apz=g*az-g1*dotprd*vz/v2
      a=apx*apx+apy*apy+apz*apz
      a=dsqrt(a)
      apx=apx/a
      apy=apy/a
      apz=apz/a
c     Compute component of B that is normal to the shock front
      dotprd=bx*apx+by*apy+bz*apz
c     Calculate components of B field parallel and perpendicular to shock normal
      bprpx=dotprd*apx
      bprpy=dotprd*apy
      bprpz=dotprd*apz
      bprp=bprpx*bprpx+bprpy*bprpy+bprpz*bprpz
      bprp=sqrt(bprp)
      bparx=bx-bprpx
      bpary=by-bprpy
      bparz=bz-bprpz
      bpar=bparx*bparx+bpary*bpary+bparz*bparz
      bpar=sqrt(bpar)
c      write(6,9999)eta,g1,vx,vy,vz,esx,esy,esz,bx,by,bz,
c     ,  bprpx,bprpy,bprpz,bparx,bpary,bparz,bpx,bpy,bpz,
c     ,  dotprd,v2,bdx,bdy,bdz
c 9999 format('*',1p25e11.3)
      return
      end
c
c     bcalc computes magnetic field component parallel and perpendicular to shock
c     front or line of sight; follows Lyutikov et al. (2003, ApJ, 597, 998)
c
      subroutine bcalc(vx,vy,vz,v,g,sx,sy,sz,bx,by,bz,
     ,  bparx,bpary,bparz,bprpx,bprpy,bprpz,bpar,bprp)
      real*8 vx,vy,vz,sx,sy,sz,s,v2,v,g,g1,g2,gp1,spx,spy,spz,
     ,  denom,spx2,spy2,spz2
c     v is velocity in units of c, s is line-of-sight or shock front unit vector
c     b is magnetic field vector
      v2=v*v
      g2=g*g
      g1=g-1.0d0
      gp1=g+1.0d0
c     Determine unit vector of shock front or l.o.s. in plasma frame
      dotprd=vx*sx+vy*sy+vz*sz
      spx=g*sx-g1*dotprd*vx/v2
      spy=g*sy-g1*dotprd*vy/v2
      spz=g*sz-g1*dotprd*vz/v2
      s=spx*spx+spy*spy+spz*spz
      s=dsqrt(s)
      spx=spx/s
      spy=spy/s
      spz=spz/s
c     Calculate components of B field parallel and perpendicular to shock front
c     or l.o.s.
      dotprd=bx*spx+by*spy+bz*spz
      bparx=dotprd*spx
      bpary=dotprd*spy
      bparz=dotprd*spz
      bpar=bparx*bparx+bpary*bpary+bparz*bparz
      bpar=sqrt(bpar)
      bprpx=bx-bparx
      bprpy=by-bpary
      bprpz=bz-bparz
      bprp=bprpx*bprpx+bprpy*bprpy+bprpz*bprpz
      bprp=sqrt(bprp)
c      write(5,9999)sx,sy,sz,spx,spy,spz,s,spx2,spy2,spz2,
c     ,  bparx,bpary,bparz,bprpx,bprpy,bprpz,bpar,bprp
c 9999 format('**',1p18e12.4)
      return
      end
c
c     vecrot rotates unit vector a by an angle psi about unit vector c along
c     a great circle to create new vector v
c
      subroutine vecrot(ax,ay,az,cx,cy,cz,psi,vx,vy,vz)
      cs=cos(psi)
      s=sin(psi)
      cs1=1.0-cs
      vx=ax*(1.0+cs1*(cx*cx-1.0))-ay*(cz*s-cs1*cx*cy)+
     ,  az*(cy*s+cs1*cx*cz)
      vy=ax*(cz*s+cs1*cx*cy)+ay*(1.0+cs1*(cy*cy-1.0))-
     ,  az*(cx*s-cs1*cy*cz)
      vz=ax*(-cy*s+cs1*cx*cz)+ay*(cx*s+cs1*cy*cz)+
     ,  az*(1.0+cs1*(cz*cz-1.0))
      return
      end
c
      subroutine scatcs(vx,vy,vz,sx,sy,sz,x,y,z,cscat)
      real*8 vx,vy,vz,sx,sy,sz,g,g1,g2,slx,sly,slz,
     ,  smx,smy,smz
c     v is velocity in units of c, s is line-of-sight unit vector
c     (x,y,z) is position of cell relative to Mach disk, cscat is
c     cosine of scattering angle in plasma frame
      v2=vx*vx+vy*vy+vz*vz
      g2=1.0/(1.0d0-v2)
      g=dsqrt(g2)
      g1=g-1.0d0
c     Determine unit vector of l.o.s. in plasma frame (slx,sly,slz)
      slx=(1.0+g1*vx*vx/v2)*sx+g1*vx*vy*sy/v2+g1*vx*vz*sz/v2
      sly=g1*vx*vy*sx/v2+(1.0+g1*vy*vy/v2)*sy+g1*vy*vz*sz/v2
      slz=g1*vx*vz*sx/v2+g1*vy*vz*sy/v2+(1.0+g1*vz*vz/v2)*sz
      s=slx*slx+sly*sly+slz*slz
      s=sqrt(s)
      slx=slx/s
      sly=sly/s
      slz=slz/s
c     Determine unit vector from Mach disk to cell in plasma frame (smx,smy,smz)
      smx=(1.0+g1*vx*vx/v2)*x+g1*vx*vy*y/v2+g1*vx*vz*z/v2
      smy=g1*vx*vy*x/v2+(1.0+g1*vy*vy/v2)*y+g1*vy*vz*z/v2
      smz=g1*vx*vz*x/v2+g1*vy*vz*y/v2+(1.0+g1*vz*vz/v2)*z
      s=smx*smx+smy*smy+smz*smz
      s=sqrt(s)
      smx=smx/s
      smy=smy/s
      smz=smz/s
c     Calculate cosine of scattering angle in plasma frame using head-on
c      approximation (photon is scattering in direction of electron motion)
      cscat=slx*smx+sly*smy+slz*smz
c      write(5,9999)sx,sy,sz,spx,spy,spz,s,spx2,spy2,spz2,
c     ,  bparx,bpary,bparz,bprpx,bprpy,bprpz,bpar,bprp
c 9999 format('**',1p18e12.4)
      return
      end
c
c  Subroutine to calculate the Klein-Nishina cross-section; from Dermer & Menon
c     (2009, High Energy Radiation from Black Holes: Gamma Rays, Cosmic Rays,
c      and Neutrinos, Princeton U. Press), eq. 6.31
c
      subroutine xseckn(epsi,epsf,g,y,xsec)
      real*8 xsec
      z=epsf/(g*epsi*y)
      xsec=(y+1.0/y-2.0*z+z*z)/(g*epsi)
      return
      end
c  Subroutine to calculate inverse Compton emission from external sources of seed photons
      function ecdust(anuf)
      common/cparm/zred1,bfield,bperp
      common/cdust/csth1,csth2,dcsth1,dcsth2,dsang,tdust
      common/cdist/gam,edist
      common/cvel/bdx,bdy,bdz,gammad,betad
      common/cseed/dnu,di,cscat
      common/crite/iwrite
      dimension gam(44),edist(44),dnu(22),di(22)
      real*8 gam,bdx,bdy,bdz,gammad,betad,dcsth1,dcsth2,dsnth1,
     ,  dsnth2,val1,val2,vala,valb,csth1,csth2,xsecc,xsecx
      ecdust=0.0
      gran=0.0
      nhi=22
      nuhi=22
c     Set the highest frequency of seed photons for scattering
      do 2011 nn=10,nhi
      if(dnu(nn).lt.anuf.and.di(nn).gt.1.0e-30) go to 2010
      nuhi=nn-1
      go to 2012
 2010 continue
 2011 continue
 2012 continue
c     Flux will be in mJy, so set x-section as (3e26/4)sigt
      s0=25.0
      homc2=8.099e-21
      ie1=1
      ef=homc2*anuf
      ef1=ef+1.0
      if(ef1.le.gam(1))go to 1
      if(ef1.ge.gam(44).or.(gam(44)/gam(1)).lt.1.01)go to 4000
      do 10 ie=2,44
      if(gam(ie).gt.ef1)go to 11
   10 continue
      ie=44
   11 ie1=ie-1
      g1=ef1
      go to 12
c     In this version, approximate that plasma velocity is along jet axis
c       and that Doppler factor of dust torus is the mean over its solid angle
c       as viewed in the plasma frame
c      tdel=1.0/(gammad*(1.0d0-betad*dcsth1))
    1 g1=gam(1)
   12 if((gam(ie1+1)/gam(ie1)).gt.1.0002)go to 13
      ie1=ie1+1
      if(ie1.ge.43)go to 4000
      go to 12
   13 vala=s0*edist(ie1)/g1
c     Loop to integrate over electron Lorentz factors (gam)
      do 3000 ie=ie1,43
      g2=gam(ie+1)
      valb=s0*edist(ie+1)/g2
      val1=0.0
      val2=0.0
      gran1=0.0
      addit=0.0
c     Set up loop 1 to integrate over incident photon frequency anui
      anumax=amin1(anuf,dnu(nuhi))
      di1=di(1)
      id=2
      anumin=anuf/((2.0d0*g1)*(g1-ef)*
     ,  (1.0d0-cscat*dsqrt(1.0d0-1.0d0/(g1*g1))))
c      anumin=0.25*anuf/((g1*g1)*(1.0-ef/g1))
      if(anumin.ge.anumax)go to 601
      if(anumin.gt.dnu(1))go to 2
      anumin=dnu(1)
      go to 5
    2 continue
      do 3 id=2,nuhi
      if(anumin.le.dnu(id))go to 4
    3 continue
      id=nuhi
    4 continue
      a=alog10(di(id-1)/di(id))/alog10(dnu(id-1)/dnu(id))
      di1=di(id-1)*(anumin/dnu(id-1))**a
    5 ide=nuhi
      if(anumax.ge.dnu(nuhi))go to 8
      do 6 idd=id,nuhi
      if(anumax.le.dnu(idd))go to 7
    6 continue
      idd=nuhi
    7 a=alog10(di(idd-1)/di(idd))/alog10(dnu(idd-1)/dnu(idd))
      die=di(idd-1)*(anumax/dnu(idd-1))**a
      ide=idd
      go to 9
    8 die=di(nuhi)
    9 continue
      anui1=anumin
      xx=1.0-ef/g1
      rat=homc2*anui1*g1*(1.0d0-cscat*dsqrt(1.0d0-1.0d0/(g1*g1)))
      call xseckn(rat,ef,g1,xx,xsecc)
      xsecc=xsecc*(g1*homc2*anui1)
      val1=xsecc*(1.0e20/anui1)*(anuf/anui1)*di1*vala
      if(val1.lt.1.0d-40)val1=0.0d0
   25 continue
c     Loop 1 to integrate over incoming photon frequency anui for lower gam value
      do 600 nu=id,ide
      anui2=dnu(nu)
      di2=di(nu)
      if(nu.lt.ide)go to 30
      anui2=anumax
   30 rat=homc2*anui2*g1*(1.0d0-cscat*dsqrt(1.0d0-1.0d0/(g1*g1)))
      call xseckn(rat,ef,g1,xx,xsecc)
      xsecc=xsecc*(g1*homc2*anui2)
      val2=xsecc*(1.0e20/anui2)*(anuf/anui2)*di2*vala
      if(val2.lt.1.0d-40)val2=0.0d0
  525 if(val1.eq.0.0d0.or.val2.eq.0.0d0)go to 845
      test=abs((anui1-anui2)/anui1)
      if(test.lt.0.001)go to 847
      ratnu=anui2/anui1
      a=1.0+dlog10(val2/val1)/alog10(ratnu)
      if(abs(a).lt.0.01.or.abs(a).gt.5.0)go to 845
      addit=val1*(ratnu**a-1.0)*anui1/a
      go to 846
  845 addit=0.5*(val1+val2)*(anui2-anui1)
  846 gran1=gran1+addit
  847 anui1=anui2
      val1=val2
      di1=di2
  600 continue
c     End anui loop 1
  601 continue
      ratg=g2/g1
      ratgl=alog10(ratg)
      valb=s0*edist(ie+1)/g2
      val1=0.0
      val2=0.0
      gran2=0.0
c     Set up loop 2 to integrate over incident photon frequency anui
      anumax=amin1(anuf,dnu(nuhi))
      di1=di(1)
      id=2
      anumin=anuf/((2.0d0*g2)*(g2-ef)*
     ,  (1.0d0-cscat*dsqrt(1.0d0-1.0d0/(g2*g2))))
c      anumin=0.25*anuf/((g2*g2)*(1.0-ef/g2))
c      write(5,6668)g1,anuf,dnu(nuhi),anumin,anumax
c 6668 format('*g1, anuf, dnu(nuhi), anumin, anumax: ',1p5e9.2)
      if(anumin.ge.anumax)go to 1601
      if(anumin.gt.dnu(1))go to 1002
      anumin=dnu(1)
      go to 1005
 1002 continue
      do 1003 id=2,nuhi
      if(anumin.le.dnu(id))go to 1004
 1003 continue
      id=nuhi
 1004 continue
      di1=0.0
      addit=0.0
      if(di(id).lt.1.0e-25.or.di(id-1).lt.1.0e-25)go to 1005
      a=alog10(di(id)/di(id-1))/alog10(dnu(id)/dnu(id-1))
      di1=di(id)*(anumin/dnu(id))**a
 1005 ide=nuhi
      if(anumax.ge.dnu(nuhi))go to 1008
      do 1006 idd=id,nuhi
      if(anumax.le.dnu(idd))go to 1007
 1006 continue
      idd=nuhi
 1007 a=alog10(di(idd-1)/di(idd))/alog10(dnu(idd-1)/dnu(idd))
      die=di(idd-1)*(anumax/dnu(idd-1))**a
      ide=idd
      go to 1009
 1008 die=di(nuhi)
 1009 continue
      anui1=anumin
      xx=1.0-ef/g2
      rat=homc2*anui1*g2*(1.0d0-cscat*dsqrt(1.0d0-1.0d0/(g2*g2)))
      call xseckn(rat,ef,g2,xx,xsecc)
      xsecc=xsecc*(g2*homc2*anui1)
      val1=xsecc*(1.0e20/anui1)*(anuf/anui1)*di1*valb
      if(val1.lt.1.0d-40)val1=0.0d0
 1025 gran2=0.0
c     Loop 2 to integrate over incoming photon frequency anui for upper gam value
      if(id.ge.ide)go to 1601
      do 1600 nu=id,ide
      anui2=dnu(nu)
      di2=di(nu)
      if(nu.lt.ide)go to 1030
      anui2=anumax
      di2=die
 1030 rat=homc2*anui2*g2*(1.0d0-cscat*dsqrt(1.0d0-1.0d0/(g2*g2)))
      call xseckn(rat,ef,g2,xx,xsecc)
      xsecc=xsecc*(g2*homc2*anui2)
      val2=xsecc*(1.0e20/anui2)*(anuf/anui2)*di2*valb
c      if(iwrite.eq.1)write(4,6666)nuhi,cscat,anumin,
c     ,  gam(40),edist(40),anui2,g1,rat,ef,xx,rat,xsecc,val2
      if(val2.lt.1.0d-40)val2=0.0d0
      if(val1.lt.1.0e-38.or.val2.lt.1.0e-38)go to 1845
      test=abs((anui1-anui2)/anui1)
      if(test.lt.0.001)go to 1847
      ratnu=anui2/anui1
      a=1.0+dlog10(val2/val1)/alog10(ratnu)
      if(abs(a).lt.0.01.or.abs(a).gt.5.0)go to 1845
      addit=val1*(ratnu**a-1.0)*anui1/a
      go to 1846
 1845 addit=0.5*(val1+val2)*(anui2-anui1)
 1846 gran2=gran2+addit
 1847 continue
c      if(iwrite.ne.1)go to 6599
c      write(4,9050)nu,id,anuf,gran2,val1,val2,vala,valb,di1,di2,
c     , anui1,anui2,g1,g2,rat,addit
 6599 continue
      anui1=anui2
      val1=val2
      di1=di2
 1600 continue
c     End anui loop 2
 1601 continue
      if(gran1.le.1.0e-28.or.gran2.le.1.0e-28)go to 2845
      if(ratgl.lt.0.0001)go to 2845
      a=1.0+alog10(gran2/gran1)/ratgl
      if(abs(a).lt.0.01.or.abs(a).gt.5.0)go to 2845
      addit=gran1*(ratg**a-1.0)*g1/a
      go to 2846
 2845 addit=0.5*(gran1+gran2)*(g2-g1)
 2846 gran=gran+addit
c      if(iwrite.eq.1)
c     , write(4,9051)gran1,val1,val2,di1,di2,anui1,anui2,anuf,
c     , rat,addit
      g1=g2
      vala=valb
 3000 continue
c     End gam loop
      ecdust=gran*1.0e-20
c      if(iwrite.eq.1)write(4,9051)anuf,gran,rat,addit
c 9050 format('*',2i8,1p15e9.2)
c 9051 format('**',1p15e9.2)
 6666 format('In ecdust',i5,2x,1p16e12.3)
 9996 format('a ',1p12e11.3)
 9997 format('b ',1p12e11.3)
 9998 format('c ',1p12e11.3)
 9999 format('d ',1p12e11.3)
 4000 return
      end
c
c  Subroutine to calculate synchrotron self-Compton emission
c    of seed photons from other cells
      function ssc(anuf)
      common/cparm/zred1,bfield,bperp
      common/cdist/gam,edist
c      common/cvel/bdx,bdy,bdz,gammad,betad
      common/cssc/dnu,di,nuhi,cscat
      common/crite/iwrite
      dimension gam(44),edist(44),dnu(68),di(68)
      real*8 gam,addit,gran,ssca,val1,val2,xsecc,xsecx
c     Prevent underflows
      if(di(nuhi).lt.1.0e-25)nuhi=nuhi-1
      ssc=0.0
      gran=0.0
      nhi=68
      nuhi=68
c     Set the highest frequency of seed photons for scattering
      do 2011 nn=20,nhi
      if(dnu(nn).lt.anuf.and.di(nn).gt.1.0e-30) go to 2010
      nuhi=nn-1
      go to 2012
 2010 continue
 2011 continue
 2012 continue
c     Flux will be in mJy, so set x-section as (3e26/8)sigt
      s0=25.0
      homc2=8.099e-21
      ie1=1
      ef=homc2*anuf
      ef1=ef+1.0
      if(ef1.le.gam(1))go to 1
      if(ef1.ge.gam(44).or.(gam(44)/gam(1)).lt.1.01)go to 4000
      do 10 ie=2,44
      if(gam(ie).gt.ef1)go to 11
   10 continue
      ie=44
   11 ie1=ie-1
      g1=ef1
      go to 12
    1 g1=gam(1)
   12 if((gam(ie1+1)/gam(ie1)).gt.1.0002)go to 13
      ie1=ie1+1
      if(ie1.ge.43)go to 4000
      go to 12
   13 vala=s0*edist(ie1)/g1
c     Loop to integrate over electron Lorentz factors (gam)
      do 3000 ie=ie1,43
      g2=gam(ie+1)
      valb=s0*edist(ie+1)/g2
      val1=0.0
      val2=0.0
      gran1=0.0
      addit=0.0
c     Set up loop 1 to integrate over incident photon frequency anui
      anumax=amin1(anuf,dnu(nuhi))
      di1=di(1)
      id=2
      anumin=anuf/((2.0d0*g1)*(g1-ef)*
     ,  (1.0d0-cscat*dsqrt(1.0d0-1.0d0/(g1*g1))))
c      anumin=0.25*anuf/((g1*g1)*(1.0-ef/g1))
      if(anumin.ge.anumax)go to 601
      if(anumin.gt.dnu(1))go to 2
      anumin=dnu(1)
      go to 5
    2 continue
      do 3 id=2,nuhi
      if(anumin.le.dnu(id))go to 4
    3 continue
      id=nuhi
    4 continue
      di1=0.0
      if(di(id).lt.1.0e-25.or.di(id-1).lt.1.0e-25)go to 5
      a=alog10(di(id)/di(id-1))/alog10(dnu(id)/dnu(id-1))
      di1=di(id)*(anumin/dnu(id))**a
    5 ide=nuhi
      if(anumax.ge.dnu(nuhi))go to 8
      do 6 idd=id,nuhi
      if(anumax.le.dnu(idd))go to 7
    6 continue
      idd=nuhi
    7 a=alog10(di(idd-1)/di(idd))/alog10(dnu(idd-1)/dnu(idd))
      die=di(idd-1)*(anumax/dnu(idd-1))**a
      ide=idd
      go to 9
    8 die=di(nuhi)
    9 continue
      anui1=anumin
      xx=1.0-ef/g1
      rat=homc2*anui1*g1*(1.0d0-cscat*dsqrt(1.0d0-1.0d0/(g1*g1)))
      call xseckn(rat,ef,g1,xx,xsecc)
      xsecc=xsecc*(g1*homc2*anui1)
      val1=xsecc*(1.0e20/anui1)*(anuf/anui1)*di1*vala
      if(val1.lt.1.0d-40)val1=0.0d0
   25 continue
c     Loop 1 to integrate over incoming photon frequency anui for lower gam value
      do 600 nu=id,ide
      anui2=dnu(nu)
      di2=di(nu)
      if(nu.lt.ide)go to 30
      anui2=anumax
   30 rat=homc2*anui2*g1*(1.0d0-cscat*dsqrt(1.0d0-1.0d0/(g1*g1)))
      call xseckn(rat,ef,g1,xx,xsecc)
      xsecc=xsecc*(g1*homc2*anui2)
      val2=xsecc*(1.0e20/anui2)*(anuf/anui2)*di2*vala
      if(val2.lt.1.0d-40)val2=0.0d0
c      if(iwrite.eq.1)write(5,6666)nuhi,cscat,anuf,anumin,
c     ,  gam(43),edist(43),anui2,g1,
c     ,  ratt,ratr,xsecx,xx,rat,xsecc,val2
  525 if(val1.eq.0.0d0.or.val2.eq.0.0d0)go to 845
      test=abs((anui1-anui2)/anui1)
      if(test.lt.0.001)go to 847
      ratnu=anui2/anui1
      a=1.0+dlog10(val2/val1)/alog10(ratnu)
      if(abs(a).lt.0.01.or.abs(a).gt.5.0)go to 845
      addit=val1*(ratnu**a-1.0)*anui1/a
      go to 846
  845 addit=0.5*(val1+val2)*(anui2-anui1)
  846 gran1=gran1+addit
  847 anui1=anui2
      val1=val2
      di1=di2
  600 continue
c     End anui loop 1
  601 continue
      ratg=g2/g1
      ratgl=alog10(ratg)
      valb=s0*edist(ie+1)/g2
      val1=0.0
      val2=0.0
      gran2=0.0
c     Set up loop 2 to integrate over incident photon frequency anui
      anumax=amin1(anuf,dnu(nuhi))
      di1=di(1)
      id=2
      anumin=anuf/((2.0d0*g2)*(g2-ef)*
     ,  (1.0d0-cscat*dsqrt(1.0d0-1.0d0/(g2*g2))))
c      anumin=0.25*anuf/((g2*g2)*(1.0-ef/g2))
c      write(5,6668)g1,anuf,dnu(nuhi),anumin,anumax
c 6668 format('*g1, anuf, dnu(nuhi), anumin, anumax: ',1p5e9.2)
      if(anumin.ge.anumax)go to 1601
      if(anumin.gt.dnu(1))go to 1002
      anumin=dnu(1)
      go to 1005
 1002 continue
      do 1003 id=2,nuhi
      if(anumin.le.dnu(id))go to 1004
 1003 continue
      id=nuhi
 1004 continue
      di1=0.0
      addit=0.0
      if(di(id).lt.1.0e-25.or.di(id-1).lt.1.0e-25)go to 1005
      a=alog10(di(id)/di(id-1))/alog10(dnu(id)/dnu(id-1))
      di1=di(id)*(anumin/dnu(id))**a
 1005 ide=nuhi
      if(anumax.ge.dnu(nuhi))go to 1008
      do 1006 idd=id,nuhi
      if(anumax.le.dnu(idd))go to 1007
 1006 continue
      idd=nuhi
 1007 a=alog10(di(idd-1)/di(idd))/alog10(dnu(idd-1)/dnu(idd))
      die=di(idd-1)*(anumax/dnu(idd-1))**a
      ide=idd
      go to 1009
 1008 die=di(nuhi)
 1009 continue
      anui1=anumin
      xx=1.0-ef/g2
      rat=homc2*anui1*g2*(1.0d0-cscat*dsqrt(1.0d0-1.0d0/(g2*g2)))
      call xseckn(rat,ef,g2,xx,xsecc)
      xsecc=xsecc*(g2*homc2*anui1)
      val1=xsecc*(1.0e20/anui1)*(anuf/anui1)*di1*valb
      if(val1.lt.1.0d-40)val1=0.0d0
 1025 gran2=0.0
c     Loop 2 to integrate over incoming photon frequency anui for upper gam value
      if(id.ge.ide)go to 1601
      do 1600 nu=id,ide
      anui2=dnu(nu)
      di2=di(nu)
      if(nu.lt.ide)go to 1030
      anui2=anumax
      di2=die
 1030 rat=homc2*anui2*g2*(1.0d0-cscat*dsqrt(1.0d0-1.0d0/(g2*g2)))
      call xseckn(rat,ef,g2,xx,xsecc)
      xsecc=xsecc*(g2*homc2*anui2)
      val2=xsecc*(1.0e20/anui2)*(anuf/anui2)*di2*valb
      if(val2.lt.1.0d-40)val2=0.0d0
      if(val1.eq.0.0.or.val2.eq.0.0)go to 1845
      test=abs((anui1-anui2)/anui1)
      if(test.lt.0.001)go to 1847
      ratnu=anui2/anui1
      a=1.0+dlog10(val2/val1)/alog10(ratnu)
      if(abs(a).lt.0.01.or.abs(a).gt.5.0)go to 1845
      addit=val1*(ratnu**a-1.0)*anui1/a
      go to 1846
 1845 addit=0.5*(val1+val2)*(anui2-anui1)
 1846 gran2=gran2+addit
 1847 continue
c      if(anuf.gt.1.0e19.and.anuf.lt.2.0e19)
c      write(5,9050)nu,id,anuf,gran2,val1,val2,vala,valb,di1,di2,
c     , anui1,anui2,g1,g2,rat,addit
      anui1=anui2
      val1=val2
      di1=di2
 1600 continue
c     End anui loop 2
 1601 continue
      if(gran1.lt.1.0e-28.or.gran2.lt.1.0e-28)go to 2845
      if(ratgl.lt.0.0001)go to 2845
      a=1.0+alog10(gran2/gran1)/ratgl
      if(abs(a).lt.0.01.or.abs(a).gt.5.0)go to 2845
      addit=gran1*(ratg**a-1.0)*g1/a
      go to 2846
 2845 addit=0.5*(gran1+gran2)*(g2-g1)
 2846 gran=gran+addit
      g1=g2
      vala=valb
 3000 continue
c     End gam loop
      ssca=gran*1.0d-20
      if(ssca.gt.1.0d-30)ssc=ssca
c      if(anuf.gt.1.0e20.and.anuf.lt.2.0e20)
c     ,write(5,9050)anuf,gran,rat,addit
c 9050 format('***',2i5,2x,1p15e9.2)
 6666 format('In ssc',i5,2x,1p16e12.3)
 4000 return
      end

c     Subroutine to round to nearest integer
      function round(f)
      integer low, high
      real*4 f, lowDiff, highDiff
      low = f
      high = f + 1.0
      lowDiff = abs(f - low)
      highDiff = abs(f - high)
      round = high;
      if(lowDiff.lt.highDiff) round = low
      return
      end

c     Subroutine to parse the command line parameters and set internal variables
c     accordingly
      function parseArgs(dummy)
      common/cinput/daysToSimulate
      common/cfixedrand/fixedRandFlag, fixedRandFileOpened, fixedRandData, fixedRandCounter
      real*4 daysToSimulate
      integer fixedRandFlag
      character(len=32) :: opt
      character(len=32) :: arg
      integer :: i
      integer :: dummy
      integer :: total
      integer :: days

      total = command_argument_count()

      do i = 1, total
         call get_command_argument(i, opt)
         select case (opt)
         case ('-d')
            call get_command_argument(i+1, arg)
            read(arg, '(i5)' ) days
            daysToSimulate = real(days)
         case ('-t')
            fixedRandFlag = 1 ! Use "fixed" set of random numbers
            write(*,*) '% Test mode ON'
         end select
      end do
      
      write(*,*) '% Days to simulate:', daysToSimulate

      parseArgs = 0
      return
      end
