      dimension spsd(16384)
      real*8 R(32768)
      real*8 stp, angle
      real*8 dat(16)
      integer i
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
      !print *, 'R: ', R
      !print *, 'R[2048]', R(2048)
      !print *, 'R[4096]', R(4096)
      !print *, 'R[8192]', R(8192)
      !print *, 'R[16384]', R(16384)

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
