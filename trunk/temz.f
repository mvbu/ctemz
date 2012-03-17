c     Turbulent, Extreme Multi-zone (TEMZ) model for blazar multi-waveband
c       light curves, SEDs, and 44(synchrotron) polarization
c
c     Program that calculates the light curve and dynamic SED of a jet with turbulent
c       plasma passing through a conical standing shock
c     This version includes variations on times as short as the time for a cell to
c       cross the shock, as well as longer-term variations according to a given PSD slope
c     Includes synchrotron radiation and IC of seed photons from a dust torus
c       and from a Mach disk
      dimension itarra(3),
c     polc(100,1140,35,250),pac(100,1140,35,250),pq(100,1140,35,250),
c     , pu(100,1140,35,250),pqcum(35),pucum(35),ai2cum(35),
     , flux(100,1141,68),gammin(100,1141),
     , gammax(100,1141),xcell(1141),gmax0(100,1141),
     , phcell(1141),rcell(1141),ycell(1141),cosph(1141),
     , sinph(1141),bperp(1141),bfield(1141),n0(100,1141),
     , zcell(100,1141),betadx(1141),betady(1141),
     , egam(100,1141,44),betadz(1141),
     , nu(68),bx(100,1141),by(100,1141),bz(100,1141),
     , delta(100,1141),enofe(100,1141,44),edist(44),fsync(68),
     , ggam(44),dustnu(22),dusti(22),ididg(1141),spsd(16384),
     , gcnt(44),igcnt(44),fsynmd(68,4000),nouter(1140),
     , fsscmd(68,4000),fmdall(68),deltmd(1140),dmd(1140),
     , tlf(100),betamx(1140),betamy(1140),
     , betamz(1140),betamr(1140),gamamr(1140),bmdx(4000),
     , bmdy(4000),bmdz(4000),bmdtot(4000),tlf1(1140),
     , flsync(100,1141,68),flcomp(100,1141,68),absorb(68,4000),
     , fsync2(68),cosmr(1140),alphmd(68,4000),dustii(110,22),
     , alfmdc(68,4000),syseed(68),scseed(68),snu(68),ssseed(68),
     , flec(100,1141,68),flssc(100,1141,68),mdd(4000),useed(110),
     , phots(68),phalph(68),icelmx(1141),imax(1141),seedpk(110),
     , abexmd(68,4000)
      character*10 dumdum
      real*8 pol,pqcum,pucum,pq,pu,pmean,polc,ai2,
     ,  pcum,tanv0,cosv0,betaup,gamup,beta,sinz,cosz,phcell,
     ,  cosph,sinph,thlos,betau,opang,tanop,cosop,sinop,
     ,  zeta,tanz,betad,betadx,betady,betadz,slos,clos,
     ,  eta,tanxi,xi,betacs,psiup,gammad,bdx,bdy,bdz,
     ,  dcsth1,dcsth2,dth1,dth2,dsnth1,dsnth2,n0ave,
     ,  betamd,betarl,betamx,betamy,betamz,betamr,gamamr,
     ,  cosmd,tanmd,cosbm,bup,bprp,gmax0,ustob,tlfact,delt,
     ,  gamb,glow,gmrat,gmratl,gmratm,ggam,tloss,t1,t2,tlmin,
     ,  eterm1,eterm2,gammax,glim,t2max,phots
      real*4 n0mean,n0,nu
      integer dstart
c      common/ci/i,j
      common/cvel/bdx,bdy,bdz,gammad,betad
      common/cparm/zred1,bfld,bperpp
      common/cdust/dcsth1,dcsth2,dsnth1,dsnth2,dsang,tdust
      common/cssc/snu,ssseed,nuhi
      common/cdist/ggam,edist
      common/cseed/dustnu,dusti
      open (2,iostat=ios, err=9000, file='temzinp.txt',
     ,   status='old')
      open (3,iostat=ios, err=9000, file='temzspec.txt',
     , status='new')
      open (4,iostat=ios, err=9000, file='temzlc.txt',
     , status='new')
      open (5,iostat=ios, err=9000, file='temzcheck.txt',
     , status='new')
c     Input file format: 1st line characters, then zred, dgpc, alpha, p,
c       bavg, n0mean, rsize, gmaxmn, gmrat, gmin
      read(2,9111)dumdum, zred, dgpc, alpha, p, bave, psdslp,
     , uratio, rsize, gmaxmn, gmrat, gmin, betaup, zeta, thlos,
     , opang,tdust,dtdist,dtrad,zdist0,vmd
      close(2)
      call itime(itarra)
      ita=itarra(1)+itarra(2)+itarra(3)
      iseed = rand (ita)
      randstart=rand(iseed)
c     icells along axial direction, jcells along transverse direction
      icells=50
      nend=10
      jcells=3*nend*(nend-1)+1
      ancol=2*nend-1
      rbound=ancol*rsize
      mdmax=1000
c     An SED will be printed out every ispecs time steps
      ispecs=1
      ispec=1
c     Set up frequencies
      do 1 inu=1,68
      phots(inu)=0.0
    1 nu(inu)=10.0**(10.0+0.25*(inu-1))
c     alpha = (s-1)/2 is the underlying spectral index, where s = slope of 
c       electron E dist.
c     gmaxmn,gmaxmx are the min/max values of initial gamma_max of
c       electrons in a cell
      gmaxmx=gmrat*gmaxmn
c     gmin is the minimum value of initial gamma of electrons
c     2p is the slope of the volume vs. initial gamma-max law over all cells,
c          V = V0*gamma_max**(-2p)
      pexp=-2.0*p
      amppsd=5.0
c     Set up compilation of distribution of initial gamma_max values
      gmrf=gmaxmx/gmaxmn
      gmrfl=alog10(gmrf)/43.0
      do 2 ie=1,44
      gfl=alog10(gmaxmn)+gmrfl*(ie-1)
      gcnt(ie)=10.0**gfl
      igcnt(ie)=0
    2 continue
      yr=3.16e7
      rad=180.0/3.14159
      pi=3.14159
      pio2=pi*0.5
      sq3=1.73205
      parsec=3.086e18
      c=3.0e10
      emc2=8.186e-7
      cc2=1.29e-9
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
c     uratio is the ratio of energy density of electrons to that of mag. field
c     set the normalization of the electron energy distribution accordingly:
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
      betaus=betaup*sinz
      gamus=1.0/sqrt(1.0-betaus**2)
c     Shock compression factor in terms of upstream speed and angle of shock zeta
c     From Cawthorne & Cobb (1990)
      eta=gamup*betaup*sinz*dsqrt(8.0d0*(betaup*sinz)**2
     , -1.0d0/gamup**2)/dsqrt(1.0d0-(betaup*cosz)**2)
c     Compression ratio of Mach disk
c       Ultra-relativistic eq. of state assumed, so compression ratio is that
c       given by Hughes, Aller, & Aller (1989, ApJ, 341, 54)
      etac=dsqrt(8.0d0*gamup**4-1.70d1*gamup**2+9.0d0)/gamup
      tanxi=(tanz**2*(3.0*betaup**2-1.0d0)-(1-betaup**2))/
     ,(tanz*(tanz**2+1.0d0+2.0d0*betaup**2))
      xi=datan(tanxi)
      betad=dsqrt((1.0d0-(betaup*cosz)**2)**2+9.0d0*
     ,(betaup**2*cosz*sinz)**2)/(3.0*betaup*sinz)
      gammad=1.0d0/dsqrt(1.0d0-betad**2)
c     Length of a cylindrical cell in pc
      zsize=2.0*rsize/tanz
      volc=pi*rsize**2*zsize
c     Length and volume of cell in plasma proper frame
      zsizep=zsize/gammad
      volcp=volc/gammad
      svmd=sqrt(vmd)
      delobs=1.0d0/(gammad*(1.0d0-betad*clos))
      dstart=mdmax+(icells+100)*delobs+25
c     Time step in observer's frame in days
      dtfact=(1.0d0-betad*clos)/(betad*clos)
      dtime=1190.0*zsize*dtfact*zred1
c      itlast=20.0/dtime
      itlast=2.0/dtime
      print *, 'itlast = ', itlast
      print *, 'dtime = ', dtime
      mdrang=0.5*(1.0/dtfact+1.0)
c     Distance of shock from axis and apex of conical jet
      tanop=dtan(opang)
      cosop=dcos(opang)
      sinop=tanop*cosop
c     Next line is specific to the selected number of cells per slice
      rshock=(2*nend-1)*rsize
c     Distance of Mach disk from z value where conical shock intersects jet boundary
      zshock=rshock/tanz
c      write(3,8888)iseed,alpha
      expon=alpha+1.0
      exp1=0.5*expon-1.0
      anorm=(1.0+alpha)/(alpha+5.0/3.0)
c     Computation of IR seed photon density from dust torus as function of distance down jet
      do 333 id=1,(icells+nend)
      zdist=zdist0+(id-nend)*zsize+zshock
c     Calculate min & max angles of dust torus in plasma frame
      dphi1=asin(dtrad/sqrt(zdist**2+dtdist**2))
      dphi2=atan(dtdist/zdist)
      dth1=dphi1+dphi2
      dth2=dphi2-dphi1
      dcsth1=-(dcos(dth1)-betad)/(1.0d0-betad*dcos(dth1))
      dcsth2=-(dcos(dth2)-betad)/(1.0d0-betad*dcos(dth2))
      dsnth1=dsqrt(1.0d0-dcsth1**2)
      dsnth2=dsqrt(1.0d0-dcsth2**2)
c     Doppler factor of dust torus emission in frame of cells
      tdel=1.0/(gammad*(1.0d0-betad*dcsth1))
      dsang=2.0*pi*(dsnth2-dsnth1)
c     Calculate seed photon field from dust emission in plasma frame
c     Peak frequency of dust thermal emission for part of torus closest to shock
      seedpk(id)=5.88e10*tdust*tdel
c     Use this to set the frequency array of the seed photons
      do 3 inu=1,22
      dustnu(inu)=seedpk(id)*10**(-1.6+(inu-1)*0.1)
      dusti(inu)=seedph(dustnu(inu))
c     If too far downstream, dust torus angles are in 2nd quadrant; need to correct
      if(dusti(inu).lt.0.0)dusti(inu)=-dusti(inu)
    3 dustii(id,inu)=dusti(inu)
  333 continue
      dflux=0.0
      do 336 id=1,icells+nend
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
      zdist=zdist0+(id-nend)*zsize+zshock
      if(inu.eq.10.or.inu.eq.18)
     ,  write(5,9996)zdist,dustnu(inu),dustii(id,inu),seedpk(id)
    6 continue
      useed(id)=dflux/c
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
      write(3,6665)zred, dgpc, alpha, p, bave, uratio, rsize, gmaxmn,
     , gmrat, gmin, betaup, (zeta*rad), (thlos*rad), (opang*rad),
     , tdust, dtdist,dtrad,zdist0,useed(5),psdslp,vmd
      write(4,6665)zred, dgpc, alpha, p, bave, uratio, rsize, gmaxmn,
     , gmrat, gmin, betaup, (zeta*rad), (thlos*rad), (opang*rad),
     , tdust, dtdist,dtrad,zdist0,useed(5),psdslp,vmd
 6665 format('#redshift: ',f5.3,2x,'Distance in Gpc: ',f5.3/
     ,  '#spectral index: ',f4.2,2x,'filling factor exponent: ',f4.2/
     ,  '#mean unshocked magnetic field: ',f4.2,2x,'ratio of electron ',
     ,  'to mag. energy: ',e7.1/'#cell radius (pc): ',f6.3/
     ,  '#Min. value of gamma_max: ',f8.1,2x,'ratio of max. to min. ',
     ,  'values of gamma_max: ',f5.1/'#gamma_min: ',f6.1,2x,
     ,  'upstream velocity: ',f9.5,'c',2x,'shock angle: ',f6.3/
     ,  '#viewing angle: ',f6.3,2x,'opening angle: ',f6.3/
     ,  '#Dust temperature: ',f6.0,2x,'distance of center of dust ',
     ,  'torus from black hole: ',f3.1,' pc'/'#radius of torus: ',
     ,  f3.1,' pc   ','Distance of shock from central ',
     ,  'engine: ',f5.2,' pc'/'#Energy density of seed photons in ',
     ,  'plasma frame: ',f10.5/'#-Slope of PSD: ',f5.1,
     ,  '     Area of Mach disk relative to other cell',1pe9.2)
c     Set up variation of energy density according to PSD; see Done et al.,
c       1992, ApJ, 400, 138, eq. B1
      tinc=dtime
      !call psdsim(16384,-psdslp,-psdslp,1.0,tinc,spsd)
      call psdsim(8192,-psdslp,-psdslp,1.0,tinc,spsd)
      psdsum=0.0
      do 4998 ip=1,8192
      spexp=1.0/(0.5*expon+1.0)
      spsd(ip)=abs(spsd(ip))**(spexp)
c      spsd(ip)=1.0
      psdsum=psdsum+spsd(ip)/16384.0
 4998 continue
      do 4999 ip=1,8192
c      if(ip.lt.100)write(5,9996)spsd(ip)
c      write(5,9935)ip,spsd(ip)
c 9935 format(i8,1pe11.3)
 4999 spsd(ip)=amppsd*spsd(ip)/psdsum
      ip0=rand(0)*5000
      it=0
 1000 continue
c     Set parameters of each cell at initial time
c     There are icells rows of cells, with jcells cells per row
c     The central cell is a Mach disk with low velocity; its observed radiation
c      is ignored, but it is an  important source of IC seed photons
      write(3,6666)
      write(4,6667)
      i=1
      j=0
      do 999 nrow=-(nend-1),(nend-1)
      ncol=2*nend-(abs(nrow)+1)
      neven=mod(ncol,2)
      ncol=ncol/2
      do 998 ncell=-ncol,ncol
      if(ncell.eq.0.and.nrow.eq.0)go to 7
      if(neven.eq.0.and.ncell.eq.0)go to 997
      j=j+1
      jold=j
      xcell(j)=2.0*ncell*rsize
      ycell(j)=nrow*sq3*rsize
      go to 8
    7 j=jcells
      xcell(j)=0.0
      ycell(j)=0.0
      rcell(j)=0.0
      cosph(j)=1.0
      sinph(j)=0.0
      j=jold
      go to 797
    8 rcell(j)=sqrt(xcell(j)**2+ycell(j)**2)
      zcol=rcell(j)/tanz
      imax(j)=2.0*zcol/zsize
      if(imax(j).lt.2)imax(j)=2
c     nouter(j) = approx. no. of cells between cell of interest and observer
      nouter(j)=imax(j)
      cosph(j)=xcell(j)/rcell(j)
      sinph(j)=ycell(j)/rcell(j)
  797 continue
c      write(5,9876)j,nrow,ncell,xcell(j),ycell(j),rcell(j),
c     ,  cosph(j),sinph(j)
  997 continue
  998 continue
  999 continue
c 9876 format('j,nrow,ncell,x,y,r,cosph,sinph ',3i5,5f12.5)
      zrf=zshock
      xrf=0.0
c      zrf=zshock+zsize*(icells-1)-rcell(30)/tanz
c      xrf=xcell(30)
c
c     *** Set up Mach disk emission for earlier times ***
c
c     Compute time delay (no. of time steps) between when plasma passes Mach disk
c       and when it passes reference point of conical shock at column 30, in plasma frame
      zmd=zrf-zshock
      idelmd=zmd/(zsize/betad)
      write(5,9595)mdmax,dstart,idelmd,ip0,zshock,rshock,zrf,zmd,
     , zsize
 9595 format('mdmax = ',i7,' dstart = ',i7,' idelmd = ',i7,
     ,  '  ip0 = ',i7/'zshock, rshock, zrf, zmd, zsize ',
     ,     1p5e10.2)
      betadx(jcells)=0.0
      betady(jcells)=0.0
      betadz(jcells)=1.0d0/3.0d0
      betamd=betadz(jcells)
      gammd=1.0d0/dsqrt(1.0d0-betamd**2)
      dopref=gammad/gammd
c     Determine B vector of MD assuming random magnetic field orientation
      i=1
      do 130 md=1,(mdmax-1)
      xrand=rand(0)
      phi=2.0*pi*xrand
      xrand=rand(0)
      costh=2.0*(xrand-0.5)
      xrand=rand(0)
      sign=xrand-0.5
      sign=sign/abs(sign)
      thetab=sign*acos(costh)
      sinthb=sin(thetab)
      costhb=cos(thetab)
      sinphb=sin(phi)
      cosphb=cos(phi)
c     Compute B field components downstream of Mach disk shock
      idelay=dstart-mdmax+md-idelmd
      if(idelay.lt.1.or.idelay.gt.(8192-ip0))
     ,  write(5,9222)idelay,i,j,md,ip0,zshock,zcell(i,j)
      bavg=bave*sqrt(spsd(idelay+ip0))
      n0mean=n0ave*spsd(idelay+ip0)
      n0prev=n0ave*spsd(idelay+ip0-1)
ccccccccccc Test
c      if(idelay.gt.(dstart+40).and.idelay.lt.(dstart+51))
c     ,  n0mean=10.0*n0mean
c      if(idelay.gt.(dstart+41).and.idelay.lt.(dstart+52))
c     ,  n0prev=10.0*n0prev
cc      write(5,9333)md,idelay,dstart,idelmd,bavg,n0mean
 9333 format(4i7,1p8e10.2)
      j=jcells
      n0(i,j)=etac*n0mean
      n0prev=etac*n0prev
      bx(i,j)=bavg*(gamup/gammad)*sinthb*cosphb*(etac-
     ,  (etac-1.0)*(cosz*cosph(j))**2)
      by(i,j)=bavg*(gamup/gammad)*sinthb*sinphb*(etac-
     ,  (etac-1.0)*(cosz*sinph(j))**2)
      bz(i,j)=bavg*costhb*(etac*cosz**2+sinz**2)
      bfield(j)=sqrt(bx(i,j)**2+by(i,j)**2+bz(i,j)**2)
      bfld=bfield(j)
      bmdx(md)=bx(i,j)
      bmdy(md)=by(i,j)
      bmdz(md)=bz(i,j)
      bmdtot(md)=bfld
c     Calculate the initial maximum electron energy in the Mach disk
      bup=bavg*dsqrt(costhb**2+(gamup*sinthb)**2)
      bprp=bavg*dsqrt((costhb*sinz**2)**2+(gamup*sinthb*cosz**2)**2*
     ,  ((cosphb*cosph(j)**2)**2+(sinphb*sinph(j)**2)**2))
c     Next line relates this to direction of B field relative to shock
c      gmax0(i,j)=gmaxmn*(bprp/bup)**(2.0*pexp)
c     Next 3 lines assume a power-law distribution of gmax0, unrelated
c         to direction of B field
c     See http://mathworld.wolfram.com/RandomNumber.html for choosing
c      random numbers from a power-law distribution
c      xrand=rand(0)
c      if(pexp.eq.0.0)gmax0(i,j)=gmaxmx
c      if(pexp.lt.0.0)gmax0(i,j)=((gmaxmx**pexp-gmaxmn**pexp)*xrand+
c     ,  gmaxmn**pexp)**(1.0/pexp)
      gmax0(i,j)=0.2*gmaxmn
      gminmd=0.15*gmaxmn
c     Calculate energy distribution in the Mach disk cell
c     Compute energy density of photons in Mach disk, time delayed by 1 step
c       Ignore photons from dust torus since beaming is small in Mach disk
      mdi=md-1
      if(mdi.le.0)mdi=1
c     Value of SSC photon density for fast cooling case, from Sari and Esen (2001)
      uphmd=bmdtot(mdi)**2/(8.0*pi)*(sqrt(1.0+4.0*uratio)-1.0)/2.0
      ustob=8.0*pi*uphmd/(bfield(j))**2
      tlfact=cc2*bfield(j)**2*(1.0+ustob)
c     iend is the last slice of cells with energetic electrons
      iend=i
      delt=zsize*yr*3.26/(gammd*betamd)
      gamb=gmax0(i,j)/(1.0+tlfact*gmax0(i,j)*delt)
      if(gamb.lt.1.0)gamb=1.0
      glow=gminmd/(1.0+tlfact*gminmd*delt)
      if(glow.lt.1.0)glow=1.0
c     Use 0.99 instead of 1 to avoid singularity of N(gamma) at gmax0
      gmrat=0.99*gmax0(i,j)/gminmd
      gmratl=dlog10(gmrat)/32.0
      gmratm=dlog10(gminmd/glow)/11.0
      ibreak=0
      if(gamb.le.glow)ibreak=1
      write(5,9996)bfield(j),n0(i,j),gmax0(i,j),uphmd,ustob,tlfact,
     ,  bmdtot(mdi),delt,gamb,glow
      do 124 ie=1,44
      egam(i,j,ie)=gminmd*10.0**(gmratl*(ie-12))
      if(ie.lt.12)egam(i,j,ie)=glow*10.0**(gmratm*(ie-1))
      if(ie.eq.1)egam(i,j,ie)=glow+0.2*(glow*10.0**(gmratm)-glow)
      if(ibreak.eq.1.or.egam(i,j,ie).lt.gamb)go to 123
      egam(i,j,ie)=gamb
      ibreak=1
  123 ggam(ie)=egam(i,j,ie)
      tloss=(gmax0(i,j)-ggam(ie))/(tlfact*ggam(ie)*gmax0(i,j))
      t1=(gminmd-ggam(ie))/(tlfact*ggam(ie)*gminmd)
      t2=dmin1(tloss,delt)
      tlmin=delt-(gminmd-ggam(ie))/(tlfact*gminmd*ggam(ie))
      enofe(i,j,ie)=0.0
      eterm1=1.0-ggam(ie)*t1*tlfact
      eterm2=1.0-ggam(ie)*t2*tlfact
      if(eterm1.ge.0.0.and.eterm2.ge.0.0)
     ,enofe(i,j,ie)=n0(i,j)/((sen-1.0)*tlfact*ggam(ie)**(sen+1.0))*
     ,  (1.0d0-eterm2**(sen-1.0))
      if(ie.lt.12)enofe(i,j,ie)=n0(i,j)/((sen-1.0)*tlfact*
     ,  ggam(ie)**(sen+1.0))*(eterm1**(sen-1.0)-
     ,  eterm2**(sen-1.0))
c     Divide by delt since integral is over time
      enofe(i,j,ie)=enofe(i,j,ie)/delt
      edist(ie)=enofe(i,j,ie)
  124 continue
      bperpp=bfld*(2.0/3.0)
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
      if(fsync(inu).gt.0.0)nuhi=inu
      if(sstau.lt.0.01)go to 126
      srcfn=fsync(inu)/sstau
      if(sstau.gt.5.0)fsync(inu)=srcfn
      if(sstau.gt.0.01.and.sstau.le.5.0)fsync(inu)=
     ,  srcfn*(1.0-exp(-sstau))
  126 continue
c     Now calculate synchrotron emission seen by other cells
c     Need to add a lower frequency to get spectral index of nu(1)
      if(inu.gt.1)go to 127
      fq1=0.98*nu(1)
      fsyn1=ajnu(fq1/dopref)*dopref**2*bperpp
      absrb1=1.02e4*(sen+2.0)*akapnu(fq1/dopref)*bperpp/(fq1/dopref)**2
  127 continue
      fsynmd(inu,md)=ajnu(restnu/dopref)*dopref**2*bperpp
      absorb(inu,md)=1.02e4*(sen+2.0)*akapnu(restnu/dopref)*
     ,   bperpp/nu(inu)**2
      alphmd(inu,md)=10.0
      abexmd(inu,md)=1.7
      if(fsynmd(inu,md).gt.0.0.and.fsyn1.gt.0.0)
     ,alphmd(inu,md)=-alog10(fsynmd(inu,md)/fsyn1)/alog10(restnu/fq1)
      if(absorb(inu,md).gt.0.0.and.absrb1.gt.0.0)
     ,abexmd(inu,md)=-alog10(fsynmd(inu,md)/fsyn1)/alog10(restnu/fq1)
c      write(5,9994)md,inu,restnu,fsync(inu),(restnu/dopref),
c     ,   fsynmd(inu,md),sstau,
c     ,   absorb(inu,md),alphmd(inu,md),bperpp,dopref
      fq1=restnu
      fsyn1=fsynmd(inu,md)
      ssseed(inu)=fsync(inu)
cc      write(5,9333)md,it,i,j,restnu,ssseed(inu),
cc     ,  fsynmd(inu,md),alphmd(inu,md)
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
      fssc1=ssc(fq1/dopref)*dopref**2
      do 129 inu=7,68
      restnu=nu(inu)
      fsscmd(inu,md)=ssc(restnu/dopref)*dopref**2
      alfmdc(inu,md)=10.0
      if(fsscmd(inu,md).gt.0.0.and.fssc1.gt.0.0)
     ,alfmdc(inu,md)=-alog10(fsscmd(inu,md)/fssc1)/alog10(restnu/fq1)
c      write(5,9333)md,it,i,j,restnu,ssseed(inu),fsscmd(inu,md),
c     ,  alfmdc(inu,md),bperpp,dopref
      fq1=restnu
      fssc1=fsscmd(inu,md)
  129 continue
  130 continue
c
c     *** End Mach disk set-up ***
c
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
      idelay=dstart+it-1+((zrf-zcell(i,j))
     ,   +(xrf-xcell(j))*slos/(betad*clos))/zsize
      if(idelay.lt.1.or.idelay.gt.(16384-ip0))
     ,  write(5,9222)idelay,i,j,md,ip0,zshock,zcell(i,j)
      bavg=bave*sqrt(spsd(idelay+ip0))
      n0mean=n0ave*spsd(idelay+ip0)
c      if(idelay.gt.(dstart+40).and.idelay.lt.(dstart+51))
c     ,  n0mean=10.0*n0mean
c      write(5,9333)it,idelay,i,j,zlos,bavg,n0mean
      phcell(j)=datan2(sinph(j),cosph(j))
c     Determine velocity vector downstream of shock
      betadx(j)=betad*cosph(j)*dsin(opang-xi)
      betady(j)=betad*sinph(j)*dsin(opang-xi)
      betadz(j)=betad*dcos(opang-xi)
c     Determine B vector of cell assuming random magnetic field orientation
   11 continue
      xrand=rand(0)
      phi=2.0*pi*xrand
      xrand=rand(0)
      costh=2.0*(xrand-0.5)
      iseed=0
      xrand=rand(0)
      sign=xrand-0.5
      sign=sign/abs(sign)
      thetab=sign*acos(costh)
      sinthb=sin(thetab)
      costhb=cos(thetab)
      sinphb=sin(phi)
      cosphb=cos(phi)
c     Compute B field components downstream of shock in the plasma frame
c       From equation that I derived by transforming from the upstream plasma
c       frame to the stationary shock frame, compressing the component of B
c       parallel to the shock, and transforming the result to the rest frame
c       of the downstream plama
      if(j.eq.jcells)go to 12
      bx(i,j)=bavg*(gamup/gammad)*sinthb*cosphb*(eta-
     ,  (eta-1.0)*(cosz*cosph(j))**2)
      by(i,j)=bavg*(gamup/gammad)*sinthb*sinphb*(eta-
     ,  (eta-1.0)*(cosz*sinph(j))**2)
      bz(i,j)=bavg*costhb*(eta*cosz**2+sinz**2)
      go to 13
c     Set field of plasma in central cell, which is a Mach disk
   12 n0(i,j)=etac*n0mean
      bx(i,j)=bavg*(gamup/gammad)*sinthb*cosphb*(etac-
     ,  (etac-1.0)*(cosz*cosph(j))**2)
      by(i,j)=bavg*(gamup/gammad)*sinthb*sinphb*(etac-
     ,  (etac-1.0)*(cosz*sinph(j))**2)
      bz(i,j)=bavg*costhb*(etac*cosz**2+sinz**2)
      bfield(j)=sqrt(bx(i,j)**2+by(i,j)**2+bz(i,j)**2)
      bfld=bfield(j)
      bmdx(md)=bx(i,j)
      bmdy(md)=by(i,j)
      bmdz(md)=bz(i,j)
      bmdtot(md)=bfld
c     Calculate the initial maximum electron energy in the Mach disk
      bup=bavg*dsqrt(costhb**2+(gamup*sinthb)**2)
      bprp=bavg*dsqrt((costhb*sinz**2)**2+(gamup*sinthb*cosz**2)**2*
     ,  ((cosphb*cosph(j)**2)**2+(sinphb*sinph(j)**2)**2))
c     Next line relates this to direction of B field relative to shock
c      gmax0(i,j)=gmaxmn*(bprp/bup)**(2.0*pexp)
c     Next 3 lines assume a power-law distribution of gmax0, unrelated
c         to direction of B field
c     See http://mathworld.wolfram.com/RandomNumber.html for choosing
c      random numbers from a power-law distribution
      xrand=rand(0)
c      if(pexp.eq.0.0)gmax0(i,j)=gmaxmx
c      if(pexp.lt.0.0)gmax0(i,j)=((gmaxmx**pexp-gmaxmn**pexp)*xrand+
c     ,  gmaxmn**pexp)**(1.0/pexp)
      gmax0(i,j)=0.2*gmaxmn
      gminmd=0.15*gmaxmn
c     Calculate energy distribution in the Mach disk cell
c     Compute energy density of photons in Mach disk, time delayed by 1 step
c       Ignore photons from dust torus since beaming is small in Mach disk
      mdi=md-1
      if(mdi.le.0)mdi=1
      uphmd=bmdtot(mdi)**2/(8.0*pi)*(sqrt(1.0+4.0*uratio)-1.0)/2.0
      ustob=8.0*pi*uphmd/(bfield(j))**2
      tlfact=cc2*bfield(j)**2*(1.0+ustob)
c     iend is the last slice of cells with energetic electrons
      iend=i
      delt=zsize*yr*3.26/(gammd*betamd)
      gamb=gmax0(i,j)/(1.0+tlfact*gmax0(i,j)*delt)
      if(gamb.lt.1.0)gamb=1.0
      glow=gminmd/(1.0+tlfact*gminmd*delt)
      if(glow.lt.1.0)glow=1.0
      gmrat=0.99*gmax0(i,j)/gminmd
      gmratl=dlog10(gmrat)/32.0
      gmratm=dlog10(gminmd/glow)/11.0
      ibreak=0
      if(gamb.le.glow)ibreak=1
cc      write(5,9996)bfield(j),n0(i,j),gmax0(i,j),useed,usdmd,ustob,tlfact,
cc     ,  delt,gamb,glow
      do 1124 ie=1,44
      egam(i,j,ie)=gminmd*10.0**(gmratl*(ie-12))
      if(ie.lt.12)egam(i,j,ie)=glow*10.0**(gmratm*(ie-1))
      if(ie.eq.1)egam(i,j,ie)=glow+0.2*(glow*10.0**(gmratm)-glow)
      if(ibreak.eq.1.or.egam(i,j,ie).lt.gamb)go to 1123
      egam(i,j,ie)=gamb
      ibreak=1
 1123 ggam(ie)=egam(i,j,ie)
      tloss=(gmax0(i,j)-ggam(ie))/(tlfact*ggam(ie)*gmax0(i,j))
      t1=(gminmd-ggam(ie))/(tlfact*ggam(ie)*gminmd)
      t2=dmin1(tloss,delt)
      tlmin=delt-(gminmd-ggam(ie))/(tlfact*gminmd*ggam(ie))
      enofe(i,j,ie)=0.0
      eterm1=1.0-ggam(ie)*t1*tlfact
      eterm2=1.0-ggam(ie)*t2*tlfact
      if(eterm1.ge.0.0.and.eterm2.ge.0.0)
     ,enofe(i,j,ie)=n0(i,j)/((sen-1.0)*tlfact*ggam(ie)**(sen+1.0))*
     ,  (1.0d0-eterm2**(sen-1.0))
      if(ie.lt.12)enofe(i,j,ie)=n0(i,j)/((sen-1.0)*tlfact*
     ,  ggam(ie)**(sen+1.0))*(eterm1**(sen-1.0)-
     ,  eterm2**(sen-1.0))
c     Divide by delt since integral is over time
      enofe(i,j,ie)=enofe(i,j,ie)/delt
      edist(ie)=enofe(i,j,ie)
 1124 continue
      bperpp=bfld*(2.0/3.0)
      nuhi=1
      do 1125 inu=1,40
      alfmdc(inu,md)=10.0
      fsscmd(inu,md)=0.0
      snu(inu)=nu(inu)
      restnu=nu(inu)
c     Synchrotron mean intensity for SSC calculation inside Mach disk
      ssabs=1.02e4*(sen+2.0)*akapnu(restnu)*bperpp/nu(inu)**2
      sstau=ssabs*parsec*rsize*svmd
      fsync(inu)=ajnu(restnu)*bperpp*rsize*svmd*(emfact*parsec)
      if(fsync(inu).gt.0.0)nuhi=inu
      if(sstau.lt.0.01)go to 1126
      srcfn=fsync(inu)/sstau
      if(sstau.gt.5.0)fsync(inu)=srcfn
      if(sstau.gt.0.01.and.sstau.le.5.0)fsync(inu)=
     ,  srcfn*(1.0-exp(-sstau))
 1126 continue
c     Now calculate synchrotron emission seen by other cells
c     Need to add a lower frequency to get spectral index of nu(1)
      if(inu.gt.1)go to 1127
      fq1=0.98*nu(1)
      fsyn1=ajnu(fq1/dopref)*dopref**2*bperpp
      absrb1=1.02e4*(sen+2.0)*akapnu(fq1/dopref)*bperpp/(fq1/dopref)**2
 1127 continue
      fsynmd(inu,md)=ajnu(restnu/dopref)*dopref**2*bperpp
      absorb(inu,md)=1.02e4*(sen+2.0)*akapnu(restnu/dopref)*
     ,   bperpp/nu(inu)**2
      alphmd(inu,md)=10.0
      abexmd(inu,md)=1.7
      if(fsynmd(inu,md).gt.0.0.and.fsyn1.gt.0.0)
     ,alphmd(inu,md)=-alog10(fsynmd(inu,md)/fsyn1)/alog10(restnu/fq1)
      if(absorb(inu,md).gt.0.0.and.absrb1.gt.0.0)
     ,abexmd(inu,md)=-alog10(fsynmd(inu,md)/fsyn1)/alog10(restnu/fq1)
c      write(5,9996)restnu,fsync(inu),fsynmd(inu,md),sstau,
c     ,   absorb(inu,md),alphmd(inu,md),bperpp
      fq1=restnu
      fsyn1=fsynmd(inu,md)
      ssseed(inu)=fsync(inu)
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
      fssc1=ssc(fq1/dopref)*dopref**2
      do 1129 inu=7,68
      restnu=nu(inu)
      fsscmd(inu,md)=ssc(restnu/dopref)*dopref**2
      alfmdc(inu,md)=10.0
      if(fsscmd(inu,md).gt.0.0.and.fssc1.gt.0.0)
     ,alfmdc(inu,md)=-alog10(fsscmd(inu,md)/fssc1)/alog10(restnu/fq1)
c      write(5,9996)restnu,fsscmd(inu,md),alfmdc(inu,md)
      fq1=restnu
      fssc1=fsscmd(inu,md)
 1129 continue
ccc   TEST
c      if(it.gt.1)go to 1131
c      do 1130 md=1,mdmax
c      write(5,9333)
c     ,  it,i,md,(idelay+md-mdmax),spsd(idelay+md-mdmax),
c     ,  fsynmd(8,md),fsynmd(20,md)
c 1130 continue
c      go to 1132
c 1131 write(5,9333)
c     ,  it,i,md,(idelay+md-mdmax),spsd(idelay+md-mdmax),
c     ,  fsynmd(8,md),fsynmd(20,md)
c 1132 continue
      go to 88
   13 bperp(j)=dsqrt(by(i,j)**2+(bz(i,j)*slos-bx(i,j)*clos)**2)
      bfield(j)=sqrt(bx(i,j)**2+by(i,j)**2+bz(i,j)**2)
      n0(i,j)=eta*n0mean
c     Calculate the initial maximum electron energy in each cell
c       from the ratio of B(perp. to shock) to B_total in shock frame
      bup=bavg*dsqrt(costhb**2+(gamup*sinthb)**2)
      bprp=bavg*dsqrt((costhb*sinz**2)**2+(gamup*sinthb*cosz**2)**2*
     ,  ((cosphb*cosph(j)**2)**2+(sinphb*sinph(j)**2)**2))
c     Next line relates this to direction of B field relative to shock
c      gmax0(i,j)=gmaxmn*(bprp/bup)**(2.0*pexp)
c     Next 3 lines assume a power-law distribution of gmax0, unrelated
c         to direction of B field
c     See http://mathworld.wolfram.com/RandomNumber.html for choosing
c      random numbers from a power-law distribution
      xrand=rand(0)
      if(pexp.eq.0.0)gmax0(i,j)=gmaxmx
      if(pexp.lt.0.0)gmax0(i,j)=((gmaxmx**pexp-gmaxmn**pexp)*xrand+
     ,  gmaxmn**pexp)**(1.0/pexp)
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
      emisco=0.0
      ecflux=0.0
      zcell(i,j)=zshock-(rcell(j)-rsize)/tanz
c     Determine Doppler factor relative to the observer
      psiup=opang
      betacs=betad*(dsin(psiup-xi)*cosph(j)*slos+dcos(psiup-xi)*clos)
      delta(i,j)=dsqrt(1.0d0-betad**2)/(1.0d0-betacs)
      bdx=betadx(j)
      bdy=betady(j)
      bdz=betadz(j)
      ecflux=0.0
      bperpp=bperp(j)
      bfld=bfield(j)
c     Determine velocity of the cell plasma relative to the MD plasma
      betamx(j)=betadx(j)/(gammd*(1.0d0-betadz(j)*betamd))
      betamy(j)=betady(j)/(gammd*(1.0d0-betadz(j)*betamd))
      betamz(j)=(betadz(j)-betamd)/(1.0d0-betadz(j)*betamd)
      betamr(j)=dsqrt(betamx(j)**2+betamy(j)**2+betamz(j)**2)
      gamamr(j)=1.0d0/dsqrt(1.0d0-betamr(j)**2)
      cosmr(j)=0.0d0
      if(zshock.eq.zcell(i,j)) go to 888
      tanmd=rcell(j)/(gamamr(j)*(zshock-zcell(i,j)))
c     Determine Doppler factor of the MD emission in the cell's frame
      cosmr(j)=1.0d0/dsqrt(tanmd**2+1.0d0)
  888 deltmd(j)=1.0d0/(gamamr(j)*(1.0d0-betamr(j)*cosmr(j)))
c     Determine time step of Mach disk seed photons from light-travel delay
      zcl=zcell(i,j)
      dmd(j)=sqrt((zshock-0.5*zsize-zcl)**2+rcell(j)**2)
      delcor=deltmd(j)/dopref
      mdmid=mdmax-(((zrf-zcl)*clos+(xrf-
     ,   xcell(j))*slos)-dmd(j))/(dtfact*zsize)
      md1=mdmid-mdrang
      if(md1.lt.1)md1=1
      md2=mdmid+mdrang
      if(md2.gt.mdmax)md2=mdmax
      if(md1.gt.md2)md1=md2
      amdrng=md2-md1+1.0
      do 3145 inu=1,68
 3145 fmdall(inu)=0.0
      do 2147 md=md1,md2
      do 3146 inu=1,68
      syseed(inu)=0.0
      scseed(inu)=0.0
 3146 ssseed(inu)=0.0
      bp2cor=0.0
c     Cosine of Mach disk's B field in frame of cell plasma, from
c       Lyutikov et al. (2003)
      term1=(bmdx(md)+bmdy(md)+bmdz(md))/(sq3*bmdtot(md))
      term2a=gamamr(j)/(bmdtot(md)*betamr(j))*
     ,  (bmdx(md)*betamx(j)+bmdy(md)*betamy(j)+bmdz(md)*betamz(j))
      term2b=(gamamr(j)/(gamamr(j)+1.0))*
     ,  (betamx(j)+betamy(j)+betamz(j))/sq3
      term3=gamamr(j)*(1.0-(betamx(j)+betamy(j)+betamz(j))/sq3)
      cosbm=(term1+term2a*(term2b-1.0))/term3
c     Correct for round-off error in case cosbm not within -1 to +1
      if(cosbm.gt.1.0)cosbm=1.0
      if(cosbm.lt.-1.0)cosbm=-1.0
      bpcorr=dsqrt(1.0d0-cosbm**2)*1.5d0
cc      write(5,9996)bpcorr,delcor,bmdtot(md),bmdx(md),bmdy(md),
cc     , bmdz(md),gamamr(j),betamx(j),betamy(j),betamz(j),cosbm
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
      sstau=ssabs*parsec*zsize
      ajofnu=fsynmd(inu,md)*bpcorr**(1.0+alphmd(inu,md))*
     ,  delcor**(2.0+alphmd(inu,md))
      ajofnu=ajofnu*zsize*(emfact*parsec)
      syseed(inu)=ajofnu
      srcfn=ajofnu/sstau
      if(sstau.gt.5.0)syseed(inu)=srcfn
      if(sstau.gt.0.1.and.sstau.le.5.0)syseed(inu)=
     ,  srcfn*(1.0-exp(-sstau))
      syseed(inu)=syseed(inu)*volc*vmd/(4.0*pi*zsize*dmd(j)**2)
      tauexp=sstau*dmd(j)/zsize
      if(tauexp.gt.15.0)syseed(inu)=0.0
      if(tauexp.le.15.0)syseed(inu)=syseed(inu)/exp(tauexp)
  145 scseed(inu)=0.0
      if(inu.lt.7)fmdall(inu)=fmdall(inu)+syseed(inu)/amdrng
c      write(5,9996)snu(inu),restnu,ssabs,sstau,ajofnu,
c     ,  fsynmd(inu,md),bpcorr,delcor,alphmd(inu,md),syseed(inu)
  146 continue
c
c     Calculate inverse Compton flux from Mach disk in frame of cell plasma
c
      do 147 inu=7,68
c     Inverse Compton mean intensity from Mach disk for 2nd-order
c       inverse Compton calculation in cells
      scseed(inu)=fsscmd(inu,md)*
     ,  delcor**(2.0+alfmdc(inu,md))
      scseed(inu)=scseed(inu)*volcp*vmd/(4.0*pi*dmd(j)**2)
      fmdall(inu)=fmdall(inu) + (syseed(inu)+scseed(inu))/amdrng
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
      if(fmdall(inu).le.0.0.or.fmdall(inu-1).le.0.0)go to 148
      aaa=alog10(fmdall(inu-1)/fmdall(inu))/alog10(snu(inu)/snu(inu-1))
      usdmd=usdmd+4.0*pi/c/(1.0-aaa)*fmdall(inu-1)*snu(inu-1)*
     ,  ((snu(inu)/snu(inu-1))**(1.0-aaa)-1.0)
  148 continue
  149 continue
  150 continue
c      if(md.gt.(mdmax+1))write(5,9232)i,j,md,md1,md2,usdmd,dmd(j),
c     ,  delcor,syseed(1),
c     ,  syseed(8),scseed(15),scseed(20),scseed(25),scseed(30)
 9232 format(2i4,3i6,1p12e10.2)
c      do 8147 inu=1,68
c      write(5,9911)j,md,nu(inu),(nu(inu)/deltmd(j)),
c     ,  fsscmd(inu,md),dmd(j),syseed(inu),scseed(inu),ssseed(inu),
c     ,  fsynmd(inu,md),alphmd(inu,md)
 9911 format('seed:  ',2i5,1p11e10.3)
c 8147 continue
c      write(5,9996)ssseed
c     Calculate energy distribution in the cell
      id=(zcell(i,j)-zshock)/zsize+nend
      if(id.lt.1)id=1
      if(id.gt.(icells+nend))write(3,9992)i,j,id,zcell(i,j),zshock,
     ,  zsize
      ustob=8.0*pi*(useed(id)+usdmd)/(bfield(j))**2
      tlfact=cc2*bfield(j)**2*(1.0+ustob)
      tlf1(j)=tlfact
      delt=zsize*yr*3.26/(gammad*betad)
      gamb=gmax0(i,j)/(1.0+tlfact*gmax0(i,j)*delt)
      glow=gmin/(1.0+tlfact*gmin*delt)
      gmrat=0.99*gmax0(i,j)/gmin
      gmratl=dlog10(gmrat)/32.0
      gmratm=dlog10(gmin/glow)/11.0
      ibreak=0
cc      write(5,9996)bfield(j),n0(i,j),gmax0(i,j),useed,usdmd,ustob,tlfact,
cc     ,  delt,gamb,glow
      do 90 ie=1,44
      egam(i,j,ie)=gmin*10.0**(gmratl*(ie-12))
      if(ie.lt.12)egam(i,j,ie)=glow*10.0**(gmratm*(ie-1))
      if(ie.eq.1)egam(i,j,ie)=glow+0.2*(glow*10.0**(gmratm)-glow)
      if(ibreak.eq.1.or.egam(i,j,ie).lt.gamb)go to 89
      egam(i,j,ie)=gamb
      ibreak=1
   89 ggam(ie)=egam(i,j,ie)
      tloss=(gmax0(i,j)-ggam(ie))/(tlfact*ggam(ie)*gmax0(i,j))
      t1=(gmin-ggam(ie))/(tlfact*ggam(ie)*gmin)
      t2=dmin1(tloss,delt)
      tlmin=delt-(gmin-ggam(ie))/(tlfact*gmin*ggam(ie))
      enofe(i,j,ie)=0.0
      eterm1=1.0-ggam(ie)*t1*tlfact
      eterm2=1.0-ggam(ie)*t2*tlfact
      if(eterm1.ge.0.0.and.eterm2.ge.0.0)
     ,enofe(i,j,ie)=n0(i,j)/((sen-1.0)*tlfact*ggam(ie)**(sen+1.0))*
     ,  (1.0d0-eterm2**(sen-1.0))
      if(ie.lt.12)enofe(i,j,ie)=n0(i,j)/((sen-1.0)*tlfact*
     ,  ggam(ie)**(sen+1.0))*(eterm1**(sen-1.0)-
     ,  eterm2**(sen-1.0))
c     Divide by delt since integral is over time
      enofe(i,j,ie)=enofe(i,j,ie)/delt
      edist(ie)=enofe(i,j,ie)
c      if(edist(ie).lt.0.0.or.edist(ie).gt.1.0e20)
c     ,  write(5,9994)i,j,ggam(ie),edist(ie),tlfact,t2,tloss,
c     ,  tlmin,t1,delt,glow,gmin,gmax0(i,j),n0(i,j)
c      write(5,9994)i,j,ggam(ie),edist(ie),tlfact,t2,tloss,
c     ,  t1,delt,glow,gmin,gmax0(i,j),n0(i,j)
   90 continue
      delt=zsize*3.26*yr/(gammad*betad)
      gammax(i,j)=gmax0(i,j)/(1.0d0+tlfact*delt*gmax0(i,j))
      gammin(i,j)=gmin/(1.0d0+tlfact*delt*gmin)
      emold=0.0
c     Start loop over observed frequencies
      ithin=0
      do 93 inu=1,22
   93 dusti(inu)=dustii(id,inu)
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
      emisco=ajnu(restnu)*bperpp*delta(i,j)**2
      fsync2(inu)=emisco*zsize*(emfact*parsec)
      if(ithin.eq.1)go to 91
      ssabs=1.02e4*(sen+2.0)*akapnu(restnu)*bperpp/
     ,   nu(inu)**2/delta(i,j)
      ssabs=ssabs*parsec*zsize*delta(i,j)
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
   91 if(inu.eq.1)go to 92
      if(emold.gt.0.0.and.fsync2(inu).gt.0.0)specin=
     ,  alog10(emold/fsync2(inu))/alog10(nu(inu)/nu(inu-1))
   92 continue
c      if(inu.eq.1)write(4,9996)gammax(i,j),
c     ,  ggam,edist,restnu,bperpp
      flsync(i,j,inu)=fsync2(inu)*(volc/zsize)*zred1/
     ,  (1.0e18*amjy*dgpc**2)*fgeom
      flux(i,j,inu)=flsync(i,j,inu)
c      if(inu.eq.19)write(5,9989)i,j,idelay,restnu,n0(i,j),bperpp,
c     ,  gammax(i,j),gammin(i,j),edist(30),edist(40)
      if(restnu.lt.1.0e14)go to 94
      spxec=0.0001
      spxssc=0.0001
      ecflux=ecdust(restnu)*3.086*volc*zred1*delta(i,j)**2/dgpc**2*
     ,   fgeom
      sscflx=ssc(restnu)*3.086*volc*zred1*delta(i,j)**2/dgpc**2*
     ,   fgeom
      taupp=0.0
      if(nu(inu).lt.1.0e22)go to 99
c     Pair production opacity calculation
c     Expression for anumin includes typical interaction angle
      anumin=(1.24e20/(nu(inu)*zred1))*1.24e20*(2.0*gammad)**2
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
      if(emeold.gt.0.0.and.ecflux.gt.0.0)spxec=
     ,  alog10(emeold/ecflux)/alog10(nu(inu)/nu(inu-1))
      if(emsold.gt.0.0.and.sscflx.gt.0.0)spxssc=
     ,  alog10(emsold/sscflx)/alog10(nu(inu)/nu(inu-1))
   94 continue
cc      if(inu.eq.53)write(5,9994)i,j,restnu,n0(i,j),gammax(i,j),
cc     ,   gammin(i,j),ecflux,sscflx,syseed(20),scseed(35)
      emold=fsync2(inu)
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
      if(icells.eq.1)go to 2200
      i=i+1
      istart=i
c     Start loop over downstream cells
      do 200 j=1,(jcells-1)
c     iend is the last slice of cells with energetic electrons
      iend=1
      ididg(j)=0
      do 200 i=istart,imax(j)
      icelmx(j)=i
      ncells=ncells+1
      if(it.gt.1)go to 110
      zcell(i,j)=zshock-(rcell(j)-rsize)/tanz
c     Determine Doppler factor relative to the observer
      psiup=opang
      betacs=betad*(dsin(psiup-xi)*cosph(j)*slos+dcos(psiup-xi)*clos)
      delta(i,j)=dsqrt(1.0d0-betad**2)/(1.0d0-betacs)
c     Determine B vector of cell assuming random magnetic field orientation
      xrand=rand(0)
      phi=2.0*pi*xrand
      xrand=rand(0)
      costh=2.0*(xrand-0.5)
      xrand=rand(0)
      sign=xrand-0.5
      sign=sign/abs(sign)
      thetab=sign*acos(costh)
      sinthb=sin(thetab)
      costhb=cos(thetab)
      sinphb=sin(phi)
      cosphb=cos(phi)
c     Compute B field components downstream of shock
      idelay=dstart+it-1-((zrf-zcell(i,j))*clos
     ,   +(xrf-xcell(j))*slos)/zsize
      if(idelay.lt.1.or.idelay.gt.(16384-ip0))
     ,  write(5,9222)idelay,i,j,md,ip0,zshock,zcell(i,j)
      bavg=bave*sqrt(spsd(idelay+ip0))
      n0mean=n0ave*spsd(idelay+ip0)
c      if(idelay.gt.(dstart+40).and.idelay.lt.(dstart+51))
c     ,  n0mean=10.0*n0mean
      n0(i,j)=eta*n0mean
      bx(i,j)=bavg*(gamup/gammad)*sinthb*cosphb*(eta-
     ,  (eta-1.0)*(cosz*cosph(j))**2)
      by(i,j)=bavg*(gamup/gammad)*sinthb*sinphb*(eta-
     ,  (eta-1.0)*(cosz*sinph(j))**2)
      bz(i,j)=bavg*costhb*(eta*cosz**2+sinz**2)
c     Calculate the initial maximum electron energy in each cell
      bup=bavg*dsqrt(costhb**2+(gamup*sinthb)**2)
      bprp=bavg*dsqrt((costhb*sinz**2)**2+(gamup*sinthb*cosz**2)**2*
     ,  ((cosphb*cosph(j)**2)**2+(sinphb*sinph(j)**2)**2))
c     Next line relates this to direction of B field relative to shock
c      gmax0(i,j)=gmaxmn*(bprp/bup)**(2.0*pexp)
c     Next 3 lines assume a power-law distribution of gmax0, unrelated
c         to direction of B field
c     See http://mathworld.wolfram.com/RandomNumber.html for choosing
c      random numbers from a power-law distribution
      xrand=rand(0)
      if(pexp.eq.0.0)gmax0(i,j)=gmaxmx
      if(pexp.lt.0.0)gmax0(i,j)=((gmaxmx**pexp-gmaxmn**pexp)*xrand+
     ,  gmaxmn**pexp)**(1.0/pexp)
      do 105 ig=1,43
  105 if(gmax0(i,j).le.gcnt(ig+1).and.gmax0(i,j).ge.gcnt(ig))
     , igcnt(ig)=igcnt(ig)+1
c     Time loop resumes here
  110 bperp(j)=dsqrt(by(i,j)**2+(bz(i,j)*slos-bx(i,j)*clos)**2)
      bfield(j)=sqrt(bx(i,j)**2+by(i,j)**2+bz(i,j)**2)
c      write(5,9996)bperp(j),bfield(j),bx(i,j),by(i,j),bz(i,j),
c     ,  slos,clos
      emisco=0.0
      ecflux=0.0
      do 111 inu=1,68
      flcomp(i,j,inu)=0.0
      flsync(i,j,inu)=0.0
      flec(i,j,inu)=0.0
      flssc(i,j,inu)=0.0
  111 flux(i,j,inu)=0.0
      bdx=betadx(j)
      bdy=betady(j)
      bdz=betadz(j)
      delt=(i-1)*zsize*yr*3.26/(gammad*betad)
      bperpp=bperp(j)
      bfld=bfield(j)
c     Determine velocity of the cell plasma relative to the MD plasma
      betamx(j)=betadx(j)/(gammd*(1.0d0-betadz(j)*betamd))
      betamy(j)=betady(j)/(gammd*(1.0d0-betadz(j)*betamd))
      betamz(j)=(betadz(j)-betamd)/(1.0d0-betadz(j)*betamd)
      betamr(j)=dsqrt(betamx(j)**2+betamy(j)**2+betamz(j)**2)
      gamamr(j)=1.0d0/dsqrt(1.0d0-betamr(j)**2)
      cosmr(j)=0.0d0
      if(zshock.eq.zcell(i,j)) go to 988
      tanmd=rcell(j)/(gamamr(j)*(zshock-zcell(i,j)))
c     Determine Doppler factor of the MD emission in the cell's frame
      cosmr(j)=1.0d0/dsqrt(tanmd**2+1.0d0)
  988 deltmd(j)=1.0d0/(gamamr(j)*(1.0d0-betamr(j)*cosmr(j)))
c     Determine time step of Mach disk seed photons from light-travel delay
      zcl=zcell(i,j)
      dmd(j)=sqrt((zshock-0.5*zsize-zcl)**2+rcell(j)**2)
      delcor=deltmd(j)/dopref
      mdmid=mdmax-(((zrf-zcl)*clos+(xrf-
     ,   xcell(j))*slos)-dmd(j))/(dtfact*zsize)
      md1=mdmid-mdrang
      if(md1.lt.1)md1=1
      md2=mdmid+mdrang
      if(md2.gt.mdmax)md2=mdmax
      if(md1.gt.md2)md1=md2
      amdrng=md2-md1+1.0
      do 4145 inu=1,68
 4145 fmdall(inu)=0.0
      do 4147 md=md1,md2
      do 4146 inu=1,68
      syseed(inu)=0.0
      scseed(inu)=0.0
 4146 ssseed(inu)=0.0
      bp2cor=0.0
c     Cosine of Mach disk's B field in frame of cell plasma, from
c       Lyutikov et al. (2003)
      term1=(bmdx(md)+bmdy(md)+bmdz(md))/(sq3*bmdtot(md))
      term2a=gamamr(j)/(bmdtot(md)*betamr(j))*
     ,  (bmdx(md)*betamx(j)+bmdy(md)*betamy(j)+bmdz(md)*betamz(j))
      term2b=(gamamr(j)/(gamamr(j)+1.0))*
     ,  (betamx(j)+betamy(j)+betamz(j))/sq3
      term3=gamamr(j)*(1.0-(betamx(j)+betamy(j)+betamz(j))/sq3)
      cosbm=(term1+term2a*(term2b-1.0))/term3
c     Correct for round-off error in case cosbm not within -1 to +1
      if(cosbm.gt.1.0)cosbm=1.0
      if(cosbm.lt.-1.0)cosbm=-1.0
      bpcorr=dsqrt(1.0d0-cosbm**2)*1.5d0
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
      ssabs=1.02e4*(sen+2.0)*akapnu(restnu)*bperpp/nu(inu)**2
      ssabs=ssabs*parsec*zsize
      ajofnu=fsynmd(inu,md)*bpcorr**(1.0+alphmd(inu,md))*
     ,  delcor**(2.0+alphmd(inu,md))
      ajofnu=ajofnu*zsize*(emfact*parsec)
      syseed(inu)=ajofnu
      srcfn=ajofnu/sstau
      if(sstau.gt.5.0)syseed(inu)=srcfn
      if(sstau.gt.0.1.and.sstau.le.5.0)syseed(inu)=
     ,  srcfn*(1.0-exp(-sstau))
      syseed(inu)=syseed(inu)*volc*vmd/(4.0*pi*zsize*dmd(j)**2)
      tauexp=sstau*dmd(j)/zsize
      if(tauexp.gt.15.0)syseed(inu)=0.0
      if(tauexp.le.15.0)syseed(inu)=syseed(inu)/exp(tauexp)
 1145 scseed(inu)=0.0
      if(inu.lt.7)fmdall(inu)=fmdall(inu)+syseed(inu)/amdrng
c      write(5,9996)snu(inu),restnu,ssabs,sstau,ajofnu,
c     ,  fsynmd(inu,md),bpcorr,delcor,alphmd(inu,md),syseed(inu)
 1146 continue
c
c     Calculate inverse Compton flux from Mach disk in frame of cell plasma
c
      do 1147 inu=7,68
c     Inverse Compton mean intensity from Mach disk for 2nd-order
c       inverse Compton calculation in cells
      scseed(inu)=fsscmd(inu,md)*
     ,  delcor**(2.0+alfmdc(inu,md))
      scseed(inu)=scseed(inu)*volcp*vmd/(4.0*pi*dmd(j)**2)
      fmdall(inu)=fmdall(inu) + (syseed(inu)+scseed(inu))/amdrng
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
      if(fmdall(inu).le.0.0.or.fmdall(inu-1).le.0.0)go to 1148
      aaa=alog10(fmdall(inu-1)/fmdall(inu))/alog10(snu(inu)/snu(inu-1))
      usdmd=usdmd+4.0*pi/c/(1.0-aaa)*fmdall(inu-1)*snu(inu-1)*
     ,  ((snu(inu)/snu(inu-1))**(1.0-aaa)-1.0)
 1148 continue
 1149 continue
 1150 continue
c      if(md.gt.(mdmax-1))write(5,9232)i,j,md,md1,md2,usdmd,dmd(j),
c     ,  delcor,syseed(1),
c     ,  syseed(8),scseed(15),scseed(20),scseed(25),scseed(30)
c     Skip flux calculation for cell if gamma_max is too low to emit at lowest frequency
c     Ratio of energy density of photons emitted by hot dust + Mach disk to
c       energy density of the magnetic field
      id=(zcell(i,j)-zshock)/zsize+nend
      if(id.lt.1)id=1
      if(id.gt.(icells+nend))write(5,9992)i,j,id,zcell(i,j),zshock,
     ,  zsize
      ustob=8.0*pi*(useed(id)+usdmd)/(bfield(j))**2
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
      delt=zsize*yr*3.26/(gammad*betad)
      glow=gammin(i-1,j)/(1.0+tlfact*gammin(i-1,j)*delt)
      gmrat=0.99*gammax(i-1,j)/gammin(i-1,j)
      gmratl=dlog10(gmrat)/32.0
      gmratm=dlog10(gammin(i-1,j)/glow)/11.0
      t2max=i*delt
      gamb=gammax(i-1,j)/(1.0+tlfact*gammax(i-1,j)*delt)
      ibreak=0
      if(gamb.le.glow.or.gamb.ge.gammax(i-1,j))ibreak=1
cc      write(5,9996)bfield(j),n0(i,j),gmax0(i,j),useed,usdmd,ustob,tlfact,
cc     ,  delt,gamb,glow
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
      eterm1=1.0-ggam(ie)*t1*tlfact
      eterm2=1.0-ggam(ie)*t2*tlfact
      if(eterm1.ge.0.0.and.eterm2.ge.0.0)
     ,enofe(i,j,ie)=n0(i,j)/((sen-1.0)*tlavg*ggam(ie)**(sen+1.0))*
     ,  (eterm1**(sen-1.0)-
     ,   eterm2**(sen-1.0))
c     Divide by cell crossing time since integral is over time
c      test1=eterm1**(sen-1.0)
c      test2=eterm2**(sen-1.0)
      enofe(i,j,ie)=enofe(i,j,ie)/delt
  188 edist(ie)=enofe(i,j,ie)
c      if(edist(ie).lt.0.0.or.edist(ie).gt.1.0e20)
c      write(5,9994)j,i,ggam(ie),edist(ie),tlfact,t2,tloss,
c     ,  tlmin,delt,glow,gmin,gamb,gmax0(i,j),n0(i,j)
c      if(i.eq.10.and.j.eq.6)
c     ,  write(5,9994)i,j,ggam(ie),edist(ie),tlfact,t2,tloss,
c     ,  tlmin,delt,glow,gmin,gamb,gmax0(i,j),n0(i,j),test1,test2
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
  193 dusti(inu)=dustii(id,inu)
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
      emisco=ajnu(restnu)*bperpp*delta(i,j)**2
      fsync2(inu)=emisco*zsize*(emfact*parsec)
      if(ithin.eq.1)go to 191
      ssabs=1.02e4*(sen+2.0)*akapnu(restnu)*bperpp/
     ,   nu(inu)**2/delta(i,j)
      ssabs=ssabs*parsec*zsize
      ssabs=ssabs*parsec*zsize
      if(ssabs.le.(0.1/ancol))ithin=1
      srcfn=fsync2(inu)/ssabs
      if(ssabs.gt.5.0)fsync2(inu)=srcfn
      if(ssabs.gt.0.1.and.ssabs.le.5.0)fsync2(inu)=
     ,  srcfn*(1.0-exp(-ssabs))
c     Absorption by other cells; zero if viewing angle >0 and cell on outer boundary
      tauexp=(nouter(j)-i)*ssabs
      if(rcell(j).gt.(0.98*rbound).and.xcell(j).le.0.0)
     ,         tauexp=0.0
      if(thlos.eq.0.0)tauexp=(nouter(j)-i)*ssabs
c      if(tauexp.gt.15.0)fsync2(inu)=0.0
c      if(tauexp.le.15.0)fsync2(inu)=fsync2(inu)/exp(tauexp)
  191 if(inu.eq.1)go to 192
      if(emold.gt.0.0.and.fsync2(inu).gt.0.0)specin=
     ,  alog10(emold/fsync2(inu))/alog10(nu(inu)/nu(inu-1))
  192 continue
      flsync(i,j,inu)=fsync2(inu)*(volc/zsize)*zred1/
     ,  (1.0e18*amjy*dgpc**2)*fgeom
      flux(i,j,inu)=flsync(i,j,inu)
      if(restnu.lt.1.0e14)go to 194
      spxec=0.0001
      spxssc=0.0001
      ecflux=ecdust(restnu)*3.086*volc*zred1*delta(i,j)**2/dgpc**2*
     ,   fgeom
      sscflx=ssc(restnu)*3.086*volc*zred1*delta(i,j)**2/dgpc**2*
     ,   fgeom
      if(nu(inu).lt.1.0e22)go to 199
      taupp=0.0
c     Pair production opacity calculation
c     Expression for anumin includes typical interaction angle
      anumin=(1.24e20/(nu(inu)*zred1))*1.24e20*(2.0*gammad)**2
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
      emold=fsync2(inu)
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
      do 300 j=1,(jcells-1)
      do 300 i=1,icelmx(j)
      tflux=tflux+flux(i,j,inu)
      tsflux=tsflux+flsync(i,j,inu)
      tcflux=tcflux+flcomp(i,j,inu)
      tecfl=tecfl+flec(i,j,inu)
      tsscfl=tsscfl+flssc(i,j,inu)
c      if(inu.eq.9)write(5,9994)j,i,flux(i,j,inu),flsync(i,j,inu),
c     ,  flcomp(i,j,inu),flec(i,j,inu),flssc(i,j,inu)
  300 continue
c      phots(inu)=tflux*amjy*
c     ,  (4*pi*dgpc**2*(zsize*1.0e18))/
c     , ((6.63e-27*nu(inu))*c*volc*jcells*zred1)
      phots(inu)=0.0
      phalph(inu)=alph+1.0
      if(tflold.eq.0.0.or.tflux.eq.0.0)go to 390
      if(inu.gt.1)alph=alog10(tflold/tflux)/
     ,alog10(nu(inu)/nu(inu-1))
      tflux=0.001*nu(inu)*tflux
      tsflux=0.001*nu(inu)*tsflux
      tcflux=0.001*nu(inu)*tcflux
      tecfl=0.001*nu(inu)*tecfl
      tsscfl=0.001*nu(inu)*tsscfl
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
      tflold=1000.0*tflux/nu(inu)
  500 continue
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
     ,  tfssc,nu(11),tfl11,ncells
c     Set up next time step by shifting physical conditions 1 slice down jet
      if(it.eq.itlast)go to 9000
      do 598 j=1,jcells
      do 598 i=istart,icelmx(j)+1
      gmax0(i,j)=gmax0(i-1,j)
      bx(i,j)=bx(i-1,j)
      by(i,j)=by(i-1,j)
      bz(i,j)=bz(i-1,j)
      n0(i,j)=n0(i-1,j)
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
 7777 format(i4,1x,i4,2f12.4,2f8.2,2x,f10.4)
 6666 format('# freq(Hz)  flux density (Jy Hz) spectral index')
 6667 format('#'/'#   time(d)     freq(Hz)  F(Jy Hz)  alpha ',
     ,   '    F(mJy)   ',
     ,   ' freq(Hz)  F(Jy Hz)  alpha     F(muJy) ',
     ,   '  freq(Hz)  F(Jy Hz)   alpha     F(nJy)    F(EC)   ',
     ,      '  F(SSC)  no. of live cells')
 8888 format(i5,f9.2)
 8889 format(/)
 9111 format(a10/f5.3/f5.3/f4.2/f4.2/f4.2/f3.1/e7.1/f5.3/f7.1/
     ,f5.1/f6.1/d9.5/d6.1/d5.1/d6.2/f6.1/f3.1/f3.1/f3.1/e7.2)
 9222 format('idelay out of bounds ',5i6,1p3e11.3)
 9989 format(i6,2x,i6,2x,i6,1p16e10.2)
 9990 format(///'Time = ',f8.2,' days'/'  freq(Hz))',2x,
     , 'Ftot(Jy Hz)  sp. index    Fsynch',4x,'     F(EC)',
     , 4x,'F(SSC-MD)')
 9991 format(i5,f8.2,2x,1p16e10.2,1x,i5)
 9992 format('** ',3i5,1p12e12.4)
 9993 format(e12.4,1x,i10,1x,2f8.3)
 9994 format(i5,2x,i5,1p16e12.4)
 9995 format('E dist',2x,1p2e11.4)
 9996 format(1p10e12.4)
 9997 format(3e10.2,2x,1p10e12.4)
 9998 format('# mean, std. dev. pol and EVPA',4f12.5/
     , '# mean, std. dev. pol., layers 1-7',2f12.5)
 9999 format(10f10.3)
      close (4, status='keep')
      close (3, status='keep')
 9000 stop
      end
c  Subroutine to calculate inverse Compton emission from external sources of seed photons
      function ecdust(anuf)
      common/cparm/zred1,bfield,bperp
      common/cdust/dcsth1,dcsth2,dsnth1,dsnth2,dsang,tdust
      common/cdist/gam,edist
      common/cvel/bdx,bdy,bdz,gammad,betad
      common/cseed/dnu,di
      dimension gam(44),edist(44),dnu(22),di(22)
      real*8 gam,bdx,bdy,bdz,gammad,betad,dcsth1,dcsth2,dsnth1,
     ,  dsnth2
      ecdust=0.0
      gran=0.0
c     Flux will be in mJy, so set x-section as (3e26/32)sigt
      s0=6.237
c     In this version, approximate that plasma velocity is along jet axis
c       and that Doppler factor of dust torus is the mean over its solid angle
c       as viewed in the plasma frame
c      tdel=1.0/(gammad*(1.0d0-betad*dcsth1))
      g1=gam(1)
      vala=s0*edist(1)/(g1*g1)
c     Loop to integration over electron Lorentz factors (gam)
      do 3000 ie=1,43
      g2=gam(ie+1)
      valb=s0*edist(ie+1)/(g2*g2)
      val1=0.0
      val2=0.0
      gran1=0.0
      addit=0.0
c     Set up loop 1 to integrate over incident photon frequency anui
      anumax=amin1(anuf,dnu(22))
      di1=di(1)
      id=2
      anumin=0.25*anuf/(g1*g1)
      if(anumin.gt.anumax)go to 601
      if(anumin.gt.dnu(1))go to 2
      anumin=dnu(1)
      go to 5
    2 continue
      do 3 id=2,22
    3 if(anumin.le.dnu(id))go to 4
      id=22
    4 continue
      a=alog10(di(id-1)/di(id))/alog10(dnu(id-1)/dnu(id))
      di1=di(id-1)*(anumin/dnu(id-1))**a
    5 ide=22
      if(anumax.ge.dnu(22))go to 8
      do 6 idd=id,22
    6 if(anumax.le.dnu(idd))go to 7
    7 a=alog10(di(idd-1)/di(idd))/alog10(dnu(idd-1)/dnu(idd))
      die=di(idd-1)*(anumax/dnu(idd-1))**a
      ide=idd
      go to 9
    8 die=di(22)
    9 continue
      anui1=anumin
      rat=anuf/(anui1*g1*g1)
      ratr=0.25*rat
   25 val1=(8.0+2.0*rat-rat*rat+4.0*rat*alog(ratr))*
     , (1.0e20/anui1)*(anuf/anui1)*di1*vala
      if(val1.lt.0.0)val1=0.0
c     Loop 1 to integrate over incoming photon frequency anui for lower gam value
      do 600 nu=id,ide
      anui2=dnu(nu)
      di2=di(nu)
      if(nu.lt.ide)go to 30
      anui2=anumax
      di2=die
   30 rat=anuf/(anui2*g1*g1)
      ratr=0.25*rat
  525 val2=(8.0+2.0*rat-rat*rat+4.0*rat*alog(ratr))*
     , (1.0e20/anui2)*(anuf/anui2)*di2*vala
      if(val2.lt.0.0)val2=0.0
      if(val1.eq.0.0.or.val2.eq.0.0)go to 845
      test=abs((anui1-anui2)/anui1)
      if(test.lt.0.001)go to 847
      ratnu=anui2/anui1
      a=1.0+alog10(val2/val1)/alog10(ratnu)
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
      valb=s0*edist(ie+1)/(g2*g2)
      val1=0.0
      val2=0.0
      gran2=0.0
c     Set up loop 2 to integrate over incident photon frequency anui
      anumax=amin1(anuf,dnu(22))
      di1=di(1)
      id=2
      anumin=0.25*anuf/(g1*g1)
      if(anumin.gt.anumax)go to 1601
      if(anumin.gt.dnu(1))go to 1002
      anumin=dnu(1)
      go to 1005
 1002 continue
      do 1003 id=2,22
 1003 if(anumin.le.dnu(id))go to 1004
      id=22
 1004 continue
      a=alog10(di(id-1)/di(id))/alog10(dnu(id-1)/dnu(id))
      di1=di(id-1)*(anumin/dnu(id-1))**a
 1005 ide=22
      if(anumax.ge.dnu(22))go to 1008
      do 1006 idd=id,22
 1006 if(anumax.le.dnu(idd))go to 1007
 1007 a=alog10(di(idd-1)/di(idd))/alog10(dnu(idd-1)/dnu(idd))
      die=di(idd-1)*(anumax/dnu(idd-1))**a
      ide=idd
      go to 1009
 1008 die=di(22)
 1009 continue
      anui1=anumin
      rat=anuf/(anui1*g2*g2)
      ratr=0.25*rat
 1025 val1=(8.0+2.0*rat-rat*rat+4.0*rat*alog(ratr))*
     , (1.0e20/anui1)*(anuf/anui1)*di1*valb
      if(val1.lt.0.0)val1=0.0
      gran2=0.0
c     Loop 2 to integrate over incoming photon frequency anui for upper gam value
      do 1600 nu=id,ide
      anui2=dnu(nu)
      di2=di(nu)
      if(nu.lt.ide)go to 1030
      anui2=anumax
      di2=die
 1030 rat=anuf/(anui2*g2*g2)
      ratr=0.25*rat
 1525 val2=(8.0d0+2.0d0*rat-rat*rat+4.0d0*rat*alog(ratr))*
     , (1.0e20/anui2)*(anuf/anui2)*di2*valb
      if(val2.lt.0.0)val2=0.0
      if(val1.eq.0.0.or.val2.eq.0.0)go to 1845
      test=abs((anui1-anui2)/anui1)
      if(test.lt.0.001)go to 1847
      ratnu=anui2/anui1
      a=1.0+alog10(val2/val1)/alog10(ratnu)
      if(abs(a).lt.0.01.or.abs(a).gt.5.0)go to 1845
      addit=val1*(ratnu**a-1.0)*anui1/a
      go to 1846
 1845 addit=0.5*(val1+val2)*(anui2-anui1)
 1846 gran2=gran2+addit
 1847 continue
c      if(anuf.gt.1.0e20.and.anuf.lt.2.0e20)
c     ,write(5,9050)anuf,gran2,val1,val2,di1,di2,anui1,anui2,
c     , rat,addit
      anui1=anui2
      val1=val2
      di1=di2
 1600 continue
c     End anui loop 2
 1601 continue
      if(gran1.eq.0.0.or.gran2.eq.0.0)go to 2845
      a=1.0+alog10(gran2/gran1)/ratgl
      if(abs(a).lt.0.01.or.abs(a).gt.5.0)go to 2845
      addit=gran1*(ratg**a-1.0)*g1/a
      go to 2846
 2845 addit=0.5*(gran1+gran2)*(g2-g1)
 2846 gran=gran+addit
c      if(anuf.gt.1.0e20)
c     , write(5,9050)gran1,val1,val2,di1,di2,anui1,anui2,anuf,
c     , rat,addit
      g1=g2
      vala=valb
 3000 continue
c     End gam loop
      ecdust=gran*1.0e-20
c      if(anuf.gt.1.0e20.and.anuf.lt.2.0e20)
c     ,write(5,9050)anuf,gran,rat,addit
 9050 format(1p15e9.2)
 4000 return
      end
c
c     Subroutine that calculates the seed photon intensity from emission by hot dust
      function seedph(f)
      common/cdust/dcsth1,dcsth2,dsnth1,dsnth2,dsang,tdust
      common/cvel/bdx,bdy,bdz,gammad,betad
      common/cfreq/freq
      real*8 bdx,bdy,bdz,gammad,betad,dcsth1,dcsth2,dsnth1,dsnth2,
     ,  sn1,sn2
      external sdgran
      freq=f
      sn1=dsnth1
      sn2=dsnth2
      call qg5(sn1,sn2,sdgran,ans)
c     Multiply by 2pi, but then divide by 4pi to get intensity
      seedph=0.5*ans
      return
      end
c
      function sdgran(sn)
      common/cdust/dcsth1,dcsth2,dsnth1,dsnth2,dsang,tdust
      common/cvel/bdx,bdy,bdz,gammad,betad
      common/cfreq/freq
      real*8 cs,sn,bdx,bdy,bdz,gammad,betad,dcsth1,dcsth2,dsnth1,dsnth2
      hok=4.80e-11
      sdgran=0.0
      cs=dsqrt(1.0d0-sn*sn)
      tdel=1.0/(gammad*(1.0d0-betad*cs))
      f=freq/tdel
      expon=hok*f/tdust
      if(expon.lt.0.01)go to 10
      if(expon.gt.5.0)go to 20
      val=(1.33e-26*f)*f/((9.0e20/f)*(exp(expon)-1.0))
      go to 25
   10 val=(2.76e-16*f)*tdust*(f/9.0e20)
      go to 25
   20 if(expon.gt.15.0)go to 26
      val=(1.33e-26*f)*f/((9.0e20/f)*exp(expon))
   25 sdgran=val*tdel*tdel
   26 return
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
      real a,b,ss,func
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
c     The remaining code is from Ritaban Chatterjee to create variations
c       according to an input PSD with slopes beta1 and beta2 at
c       variational frequencies below and above a break frequency nu_break
C**************************************************************
      subroutine psdsim(N,beta1,beta2,nu_break,t_incre1,lc_sim)
C**************************************************************
C N=number of data points in the lc, should be an integer power of 2, N<=8192
C t_incre1=increment in time in each step while resampling the simulated data. 
C This must be larger than the smallest interval between successive data points in the input light curve

      implicit none

      real*8 nu(16384),flux_s1(16384),dat(16384),R(32768)
      real*8 dataim_s1(16384),datareal_s1(16384),amp(16384)
      real*8 flux_s2(16384),dataim_s2(16384),datareal_s2(16384)
      real*8 dataim(16384),datareal(16384),flux_s(16384)
      real*8 fac_norm,fac_norm2
      real*4 ann,beta1,beta2,nu_break,t_incre1,lc_sim(16384)
      integer N,j,i,ifile,ik,nn,ISEED1,ISEED2,isign

      fac_norm=1./(N*t_incre1)
      fac_norm2=N**2./(2.*N*t_incre1)

      ISEED1=58
      ISEED2=256871
      call RNE2IN(ISEED1,ISEED2)

        call RNSTNR(R,N/2)
        call RNE2OT(ISEED1,ISEED2)
        do j=1,N/2  
            nu(j)= j*fac_norm/86400.                      
            if (nu(j) .gt. nu_break) go to 10
              flux_s(j)=nu(j)**beta1
              datareal(j)=dsqrt(fac_norm2*0.5*flux_s(j))*R(j)
            go to 20
   10       flux_s(j)=(nu(j)**beta2)
     >                      *((nu_break**beta1)/(nu_break**beta2))
              datareal(j)=dsqrt(fac_norm2*0.5*flux_s(j))*R(j)
   20 continue
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
c  Subroutine to calculate synchrotron self-Compton emission
c    of seed photons from other cells
      function ssc(anuf)
      common/cparm/zred1,bfield,bperp
      common/cdist/gam,edist
c      common/cvel/bdx,bdy,bdz,gammad,betad
      common/cssc/dnu,di,nuhi
      dimension gam(44),edist(44),dnu(68),di(68)
      real*8 gam,addit,gran,ssca
c     Prevent underflows
      if(di(nuhi).lt.1.0e-25)nuhi=nuhi-1
      ssc=0.0
      gran=0.0
c     Flux will be in mJy, so set x-section as (3e26/32)sigt
      s0=6.237
      g1=gam(1)
      vala=s0*edist(1)/(g1*g1)
c     Loop to integrate over electron Lorentz factors (gam)
c      write(5,6666)nuhi,anuf,g1,edist(1),gam(44),edist(44),
c     ,  vala,dnu(1),di(1),dnu(nuhi),di(nuhi)
 6666 format('In ssc',i5,2x,1p10e12.3)
      do 3000 ie=1,43
      g2=gam(ie+1)
      valb=s0*edist(ie+1)/(g2*g2)
      val1=0.0
      val2=0.0
      gran1=0.0
      addit=0.0
c     Set up loop 1 to integrate over incident photon frequency anui
      anumax=amin1(anuf,dnu(nuhi))
      di1=di(1)
      id=2
      anumin=0.25*anuf/(g1*g1)
c      write(5,6667)g1,anuf,dnu(nuhi),anumin,anumax
c 6667 format('g1, anuf, dnu(nuhi), anumin, anumax: ',1p5e9.2)
      if(anumin.ge.anumax)go to 601
      if(anumin.gt.dnu(1))go to 2
      anumin=dnu(1)
      go to 5
    2 continue
      do 3 id=2,nuhi
    3 if(anumin.le.dnu(id))go to 4
      id=nuhi
    4 continue
      di1=0.0
      if(di(id).lt.1.0e-25.or.di(id-1).lt.1.0e-25)go to 5
      a=alog10(di(id)/di(id-1))/alog10(dnu(id)/dnu(id-1))
      di1=di(id)*(anumin/dnu(id))**a
    5 ide=nuhi
      if(anumax.ge.dnu(nuhi))go to 8
      do 6 idd=id,nuhi
    6 if(anumax.le.dnu(idd))go to 7
    7 a=alog10(di(idd-1)/di(idd))/alog10(dnu(idd-1)/dnu(idd))
      die=di(idd-1)*(anumax/dnu(idd-1))**a
      ide=idd
      go to 9
    8 die=di(nuhi)
    9 continue
      anui1=anumin
      rat=anuf/(anui1*g1*g1)
      ratr=0.25*rat
   25 val1=(8.0+2.0*rat-rat*rat+4.0*rat*alog(ratr))*
     , (1.0e20/anui1)*(anuf/anui1)*di1*vala
      if(val1.lt.1.0e-28)val1=0.0
c     Loop 1 to integrate over incoming photon frequency anui for lower gam value
      do 600 nu=id,ide
      anui2=dnu(nu)
      di2=di(nu)
      if(nu.lt.ide)go to 30
      anui2=anumax
      di2=die
   30 rat=anuf/(anui2*g1*g1)
      ratr=0.25*rat
  525 val2=(8.0+2.0*rat-rat*rat+4.0*rat*alog(ratr))*
     , (1.0e20/anui2)*(anuf/anui2)*di2*vala
      if(val2.lt.1.0e-28)val2=0.0
      if(val1.eq.0.0.or.val2.eq.0.0)go to 845
      test=abs((anui1-anui2)/anui1)
      if(test.lt.0.001)go to 847
      ratnu=anui2/anui1
      a=1.0+alog10(val2/val1)/alog10(ratnu)
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
      valb=s0*edist(ie+1)/(g2*g2)
      val1=0.0
      val2=0.0
      gran2=0.0
c     Set up loop 2 to integrate over incident photon frequency anui
      anumax=amin1(anuf,dnu(nuhi))
      di1=di(1)
      id=2
      anumin=0.25*anuf/(g1*g1)
c      write(5,6668)g1,anuf,dnu(nuhi),anumin,anumax
c 6668 format('*g1, anuf, dnu(nuhi), anumin, anumax: ',1p5e9.2)
      if(anumin.ge.anumax)go to 1601
      if(anumin.gt.dnu(1))go to 1002
      anumin=dnu(1)
      go to 1005
 1002 continue
      do 1003 id=2,nuhi
 1003 if(anumin.le.dnu(id))go to 1004
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
 1006 if(anumax.le.dnu(idd))go to 1007
 1007 a=alog10(di(idd-1)/di(idd))/alog10(dnu(idd-1)/dnu(idd))
      die=di(idd-1)*(anumax/dnu(idd-1))**a
      ide=idd
      go to 1009
 1008 die=di(nuhi)
 1009 continue
      anui1=anumin
      rat=anuf/(anui1*g2*g2)
      ratr=0.25*rat
 1025 val1=(8.0+2.0*rat-rat*rat+4.0*rat*alog(ratr))*
     , (1.0e20/anui1)*(anuf/anui1)*di1*valb
      if(val1.lt.1.0e-28)val1=0.0
      gran2=0.0
c     Loop 2 to integrate over incoming photon frequency anui for upper gam value
      do 1600 nu=id,ide
      anui2=dnu(nu)
      di2=di(nu)
      if(nu.lt.ide)go to 1030
      anui2=anumax
      di2=die
 1030 rat=anuf/(anui2*g2*g2)
      ratr=0.25*rat
 1525 val2=(8.0d0+2.0d0*rat-rat*rat+4.0d0*rat*alog(ratr))*
     , (1.0e20/anui2)*(anuf/anui2)*di2*valb
      if(val2.lt.1.0e-28)val2=0.0
      if(val1.eq.0.0.or.val2.eq.0.0)go to 1845
      test=abs((anui1-anui2)/anui1)
      if(test.lt.0.001)go to 1847
      ratnu=anui2/anui1
      a=1.0+alog10(val2/val1)/alog10(ratnu)
      if(abs(a).lt.0.01.or.abs(a).gt.5.0)go to 1845
      addit=val1*(ratnu**a-1.0)*anui1/a
      go to 1846
 1845 addit=0.5*(val1+val2)*(anui2-anui1)
 1846 gran2=gran2+addit
 1847 continue
c      if(anuf.gt.1.0e19.and.anuf.lt.2.0e19)
c      write(5,9050)nu,id,anuf,gran2,val1,val2,vala,valb,di1,di2,
c     , anui1,anui2,g1,g2,rat,addit
      if(val2.lt.1.0e-28)go to 1601
      anui1=anui2
      val1=val2
      di1=di2
 1600 continue
c     End anui loop 2
 1601 continue
      if(gran1.lt.1.0e-28.or.gran2.lt.1.0e-28)go to 2845
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
 9050 format('***',2i5,2x,1p15e9.2)
 4000 return
      end
