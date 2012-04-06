      dimension spsd(16384)
      dimension edist(44), ggam(44)
      real*8 R(32768)
      real*8 stp, angle
      real*8 dat(16)
      real*4 anu
      real*8 ggam
      integer i
      common/cparm/zred1,bfield,b
      common/cdist/ggam,edist
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

