      dimension spsd(16384)
      dimension edist(44), ggam(44), dustnu(22), dusti(22)
      dimension snu(68), ssseed(68)
      dimension ggam_ecdust(44), edist_ecdust(44), edist_ssc(44)
      real*8 R(32768)
      real*8 stp, angle
      real*8 dat(16)
      real*4 anu
      real*8 ggam
      real*8 clos, slos
      integer i
      common/cparm/zred1,bfield,b
      common/cdist/ggam,edist
      common/cseed/dustnu, dusti
      common/cssc/snu, ssseed, nuhi
      common/cvel/bdx,bdy,bdz,gammad,betad
      real*8 bdx,bdy,bdz,gammad,betad
      data ggam  / 1.123009204864502, 1.6150462627410889, 2.6083743572235107, 4.2126455307006836, 6.803617000579834, 10.988156318664551, 17.746381759643555, 28.661226272583008, 46.289207458496094, 74.759208679199219, 120.73958587646484, 195.00001525878906, 196.69918823242188, 198.41317749023438, 200.14210510253906, 201.88607788085938, 203.645263671875, 205.41976928710938, 207.20974731445312, 209.01531982421875, 210.83662414550781, 212.67379760742188, 214.5269775390625, 216.39631652832031, 218.28193664550781, 220.18399047851562, 222.10261535644531, 224.03794860839844, 225.99015808105469, 227.95938110351562, 229.94575500488281, 231.94944763183594, 233.97059631347656, 236.00935363769531, 238.06587219238281, 240.14031982421875, 242.23283386230469, 244.34358215332031, 246.47273254394531, 248.62042236328125, 250.78683471679688, 252.97213745117188, 255.17646789550781, 257.39999389648438 /
      data edist / 207.588089, 100.368889, 38.4794807, 14.7522821, 5.65573835, 2.16830015, 0.831283987, 0.318698138, 0.122182652, 0.0468424521, 0.0179584846, 0.00688493345, 0.00650044205, 0.00613163831, 0.00577794248, 0.0054388023, 0.00511366921, 0.00480203005, 0.00450337958, 0.0042172363, 0.0039431327, 0.00368061941, 0.00342926243, 0.00318864221, 0.0029583571, 0.00273801689, 0.00252724718, 0.00232568663, 0.00213298434, 0.00194880483, 0.00177282514, 0.00160473085, 0.00144422159, 0.00129100704, 0.00114480685, 0.00100535038, 0.000872378412, 0.000745639321, 0.000624891079, 0.000509901613, 0.000400445395, 0.000296305661, 0.000197275149, 0.000103152641 /
      N = 16384
      stp = 6.28318530717959d0/dble(8-1d0)
      ISEED1=58
      ISEED2=256871
      print *, 'ISEED1, ISEED2: ', ISEED1, ISEED2
      call RNE2IN(ISEED1,ISEED2)
      print *, 'ISEED1, ISEED2: ', ISEED1, ISEED2
      call RNSTNR(R,N/2)
      call RNE2OT(ISEED1,ISEED2)
      print *, 'ISEED1, ISEED2: ', ISEED1, ISEED2

      NN=8
      do i=1, NN, 1 
         ! set complex part to zero
         dat(i*2) = 0.0d0
         ! set real part to a sine wave
         angle = dble(i-1)*stp
         dat(i*2-1) = sin(angle)
      end do

      !call printcomplex(dat, NN)
      call four1(dat, NN, 1)
      !call printcomplex(dat, NN)

      ! Now test psdsim
      psdslp = 1.7
      tinc= 1.93893981
      !call psdsim(16384,-psdslp,-psdslp,1.0,tinc,spsd)
      call psdsim(8192,-psdslp,-psdslp,1.0,tinc,spsd)
      print *, 'spsd[   1]', spsd(1)
      print *, 'spsd[2048]', spsd(2048)
      print *, 'spsd[4096]', spsd(4096)
      print *, 'spsd[8192]', spsd(8192)
      print *, 'spsd[16384]', spsd(16384)

      ! Test akapnu()
      b = 12.7397528
      anu = 1e10
      akapnuResult = akapnu(anu);
      print *, 'akapnuResult', akapnuResult
      print *, 'expected',  1.16971385

      ! Test ssc()
      data edist_ssc / 443.1716, 214.273575, 82.1483231, 31.4940643, 12.074213, 4.62901831, 1.77467537, 0.680376053, 0.260842919, 0.10000211, 0.0383388586, 0.0146983732, 0.013877538, 0.0130901923, 0.0123351021, 0.0116110835, 0.0109169716, 0.0102516655, 0.00961408857, 0.00900321081, 0.00841803849, 0.00785760954, 0.00732099777, 0.00680730632, 0.00631567976, 0.00584528456, 0.00539532024, 0.00496501662, 0.00455362396, 0.00416042656, 0.00378473452, 0.00342587661, 0.00308321184, 0.00275612017, 0.00244400324, 0.00214628293, 0.00186240615, 0.00159183599, 0.00133405521, 0.00108856894, 0.000854895043, 0.000632571289, 0.000421154924, 0.000220216491 /
      data snu / 1e+10, 1.77827942e+10, 3.16227768e+10, 5.62341315e+10, 9.9999998e+10, 1.77827938e+11, 3.16227781e+11, 5.62341347e+11, 9.99999996e+11, 1.77827938e+12, 3.16227768e+12, 5.62341347e+12, 9.99999983e+12, 1.77827941e+13, 3.16227768e+13, 5.62341326e+13, 1e+14, 1.77827937e+14, 3.16227772e+14, 5.6234131e+14, 9.99999987e+14, 1.7782794e+15, 3.16227758e+15, 5.62341303e+15, 1.00000003e+16, 1.77827945e+16, 3.16227769e+16, 5.62341314e+16, 9.99999984e+16, 1.77827933e+17, 3.16227756e+17, 5.62341339e+17, 9.99999984e+17, 1.77827946e+18, 3.16227763e+18, 5.62341319e+18, 9.99999998e+18, 1.77827941e+19, 3.16227768e+19, 5.62341352e+19, 1.00000002e+20, 1.77827939e+20, 3.16227777e+20, 5.62341334e+20, 1.00000002e+21, 1.77827939e+21, 3.16227763e+21, 5.62341348e+21, 9.99999978e+21, 1.77827941e+22, 3.16227769e+22, 5.62341337e+22, 9.99999978e+22, 1.77827946e+23, 3.16227778e+23, 5.62341319e+23, 1.00000001e+24, 1.77827939e+24, 3.16227749e+24, 5.62341326e+24, 9.99999956e+24, 1.77827939e+25, 3.16227761e+25, 5.62341338e+25, 1.00000003e+26, 1.77827941e+26, 3.16227774e+26, 5.62341319e+26 /
      data ssseed  / 1.03139553e-06, 4.12530653e-06, 1.72026703e-05, 5.10973186e-05, 5.7645284e-05, 4.3994758e-05, 2.86946997e-05, 1.75751247e-05, 9.08344646e-06, 4.14674287e-06, 1.46940556e-06, 3.13423982e-07, 2.33193926e-08, 3.34458017e-10, 4.20477265e-13, 4.55808116e-18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /

      ! data statement only works once for a given variable
      ! so have to create new arrays and copy them into the common block
      do i=1,44
         edist(i) = edist_ssc(i)
      end do

      nuhi = 16
      anuf = 2.86500065e+10
      sscResult = ssc(anuf)
      print *, 'sscResult', sscResult
      print *, 'expected',  0.024821050307809192
      
      ! Test ecdust()
      data ggam_ecdust / 0.87952703237533569, 1.6948215961456299, 3.2668538093566895, 6.297025203704834, 12.137832641601562, 23.396284103393555, 45.097515106201172, 86.927726745605469, 167.55757141113281, 322.97564697265625, 622.55181884765625, 1200, 1219.662841796875, 1239.6478271484375, 1259.9603271484375, 1280.6055908203125, 1301.5892333984375, 1322.9166259765625, 1344.5936279296875, 1366.6256103515625, 1389.018798828125, 1411.77880859375, 1434.9117431640625, 1458.423828125, 1482.321044921875, 1506.60986328125, 1531.2967529296875, 1556.3880615234375, 1581.890625, 1607.8109130859375, 1634.156005859375, 1660.932861328125, 1688.1483154296875, 1715.809814453125, 1743.924560546875, 1772.5, 1801.5435791015625, 1831.0631103515625, 1861.0662841796875, 1891.5611572265625, 1922.5556640625, 1954.05810546875, 1986.07666015625, 2018.619873046875 /
      data edist_ecdust / 54.3787384, 14.6458302, 3.94187236, 1.06094074, 0.285548359, 0.0768543035, 0.0206850581, 0.00556730898, 0.00149842107, 0.000403294631, 0.00010854529, 2.92145687e-05, 2.6896425e-05, 2.47395292e-05, 2.27332275e-05, 2.0867581e-05, 1.91332601e-05, 1.75215719e-05, 1.60243453e-05, 1.46339853e-05, 1.33433377e-05, 1.21457533e-05, 1.10349883e-05, 1.00052102e-05, 9.05097568e-06, 8.16718057e-06, 7.34905689e-06, 6.59215675e-06, 5.89230831e-06, 5.24562711e-06, 4.64846926e-06, 4.09743552e-06, 3.58935517e-06, 3.12125439e-06, 2.69036309e-06, 2.29409102e-06, 1.93002074e-06, 1.5958916e-06, 1.28959528e-06, 1.00916122e-06, 7.52751987e-07, 5.18650722e-07, 3.05256975e-07, 1.11075792e-07 /
      data dustnu / 1.4985991e+11, 1.88662448e+11, 2.37511967e+11, 2.99009901e+11, 3.76431051e+11, 4.73898648e+11, 5.96603109e+11, 7.51078736e+11, 9.45552163e+11, 1.19037965e+12, 1.49859913e+12, 1.88662468e+12, 2.3751198e+12, 2.99009861e+12, 3.76431051e+12, 4.7389868e+12, 5.96603096e+12, 7.51078788e+12, 9.45552255e+12, 1.19037955e+13, 1.49859913e+13, 1.8866252e+13 /
      data dusti / 9.51800444e-12, 1.48611574e-11, 2.31122291e-11, 3.57641555e-11, 5.49886907e-11, 8.38589476e-11, 1.26555683e-10, 1.88442997e-10, 2.75783063e-10, 3.94694805e-10, 5.4881516e-10, 7.35197514e-10, 9.38718325e-10, 1.12717435e-09, 1.25207189e-09, 1.26118194e-09, 1.12446863e-09, 8.61036353e-10, 5.45080592e-10, 2.72187495e-10, 1.00743747e-10, 2.54984783e-11 /

      ! data statement only works once for a given variable
      ! so have to create new arrays and copy them into the common block
      do i=1,44
         ggam(i) = ggam_ecdust(i)
         edist(i) = edist_ecdust(i)
      end do

      anuf = 1.70602292e+14
      ecdustResult = ecdust(anuf)
      print *, 'ecdustResult', ecdustResult
      print *, 'expected',  2.84165187e-7

      ! test polcalc()
      b = 1.05397856
      bx = 0.898704171
      by = -.54515934
      bz = -.0774780065
      clos = 0.9999862922696503
      slos = 0.0052359595870788813
      bdx = 0.021429364179171097
      bdy = 0.033405057415456803
      bdz = 0.99540335416720571
      expectedChi = 2.61113882
      call polcalc(b, bx, by, bz, clos, slos, chi)
      print *, 'polcalcResult', chi
      print *, 'expected',  expectedChi
      
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
       SUBROUTINE four1(dat,nn,isign)
C****************************************************************

       INTEGER isign,nn
       REAL*8 dat(2*nn)
c  Replaces dat(1:2*nn) by its discrete Fourier transform, 
c  if isign is input as 1; or replaces
c  dat(1:2*nn) by nn times its inverse discrete Fourier 
c  transform, if isign is input as -1.
c  dat is a complex array of length nn or, equivalently, 
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
            tempr=dat(j) ! Exchange the two complex numbers.
            tempi=dat(j+1)
            dat(j)=dat(i)
            dat(j+1)=dat(i+1)
            dat(i)=tempr
            dat(i+1)=tempi
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
                tempr=sngl(wr)*dat(j)-sngl(wi)*dat(j+1)
                tempi=sngl(wr)*dat(j+1)+sngl(wi)*dat(j)
                dat(j)=dat(i)-tempr
                dat(j+1)=dat(i+1)-tempi
                dat(i)=dat(i)+tempr
                dat(i+1)=dat(i+1)+tempi
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

C**************************************************************
      subroutine psdsim(N,beta1,beta2,nu_break,t_incre1,lc_sim)
C**************************************************************
C N=number of data points in the lc, should be an integer power of 2, N<=8192
C t_incre1=increment in time in each step while resampling the simulated data. 
C This must be larger than the smallest interval between successive data points in the input light curve

      implicit none

      real*8 nu(16384),flux_s1(16384),dat(16384),R(32768)
      real*8 dataim_s1(16384),datareal_s1(16384),amp(16384)
      real*8 fluxs2(16384),dataim_s2(16384),datareal_s2(16384)
      real*8 dataim(16384),datareal(16384),flux_s(16384)
      real*8 fac_norm,fac_norm2
      real*4 ann,beta1,beta2,nu_break,t_incre1,lc_sim(16384)
      integer N,j,i,ifile,ik,nn,ISEED1,ISEED2,isign

      fac_norm=1./(N*t_incre1)
      fac_norm2=N**2./(2.*N*t_incre1)

      ISEED1=58
      ISEED2=256871
      call RNE2IN(ISEED1,ISEED2)
        print *, '1) ISEED1 ', ISEED1, ' ISEED2 ', ISEED2

        call RNSTNR(R,N/2)
        call RNE2OT(ISEED1,ISEED2)
        print *, '2) ISEED1 ', ISEED1, ' ISEED2 ', ISEED2
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
        print *, '3) ISEED1 ', ISEED1, ' ISEED2 ', ISEED2

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

      SUBROUTINE printcomplex(dat,nn)
      INTEGER nn
      REAL*8 dat(2*nn)
      integer i

      do i=1, nn, 1
         print *, dat(2*i-1), dat(2*i)
      end do
      return
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

      function ecdust(anuf)
      common/cparm/zred1,bfield,bperp ! Not used as far as I can tell
      common/cdust/dcsth1,dcsth2,dsnth1,dsnth2,dsang,tdust ! Not used AFAICT
      common/cdist/gam,edist
      common/cvel/bdx,bdy,bdz,gammad,betad ! Not used AFAICT
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

      !
      ! ssc()
      !
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
c
c     polcalc computes the polarization angle chi based on relations
c       derived by Lyutikov et al. 2005, MNRAS, 360, 869
c
      Subroutine polcalc(b,bx,by,bz,clos,slos,chi)
      common/cvel/bdx,bdy,bdz,gammad,betad
      real*8 clos,slos,bdx,bdy,bdz,gammad,betad,term,ndq,q2
      bxh=bx/b
      byh=by/b
      bzh=bz/b
      term=byh+(bdy*bxh-bdx*byh)*slos+(bdy*bzh-bdz*byh)*clos 
      ndq=(bxh+(bdx*bzh-bdz*bxh)*clos)*slos +
     ,  (bzh+(bdz*bxh-bdx*bzh)*slos)*clos
      q2=(bxh+(bdx*bzh-bdz*bxh)*clos)**2+
     ,   (byh+(bdy*bxh-bdx*byh)*slos+(bdy*bzh-bdz*byh)*clos)**2 +
     ,   (bzh+(bdz*bxh-bdx*bzh)*slos)**2
      coschi=term/dsqrt(q2-ndq*ndq)
      chi=acos(coschi)
c      write(5,9999)b,bx,by,bz,clos,slos,bdx,bdy,bdz,
c     ,   gammad,term,ndq,q2,coschi,chi
c 9999 format(1p15e10.2)
      return
      end
