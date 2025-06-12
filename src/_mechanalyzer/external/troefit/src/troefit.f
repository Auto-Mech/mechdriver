
      Program troefit

c Program to fit
c (i) A0,n0,E0 in k0=A0 T^n0 exp(-E0/T)
c (ii) Ainf,ninf,Einf in kinf=Ainf T^ninf exp(-Einf/T)
c (iii) alpha,T*,T**,T*** in fcent = alpha*exp(-T/T*)+(1-alpha)*exp(-T/T***)+exp(-T**/T) 
c with the parameters specified in that order.
c for some reason A0 and Ainf are input as natural logs here
C also alpha is input as alpha=exp(-a0^2).
c also CHEMKIN paramters are in the order alpha, T***, T*, T**
c the chemkin manual sometimes forgets the minus sign in exp(-T/T*)

      implicit real*8 (a-h,o-z)
      parameter (ma=10)
      dimension pin(300),dkin(300),sig(300)
      dimension a0(ma),a(ma),dkda(ma)
      dimension ilist(ma)
      dimension covar(300,300),alpha(300,300)
      common /tp/ temp(300),pres(300)
      external funcs

      OPEN(UNIT=5,STATUS='OLD',FILE='troefit.dat')
      OPEN(UNIT=6,STATUS='unknown',FILE='troefit.out')

      read (5,*) converg1,converg2
      read (5,*) (ilist(ii),ii=1,ma)
      read (5,*) iinp
      if (iinp.eq.1) then
         read (5,*) (a0(ii),ii=1,ma)
      endif
      tmin = 1.0d30
      pmin = 0.0d30
      tmax = 0.0d0
      pmax = 0.0d0
      itp = 0
      read (5,*) nt
      do it = 1 , nt
         read (5,*) temper,np
c        write (6,*) 'temp',temper,np
         do ip = 1 , np
            itp = itp+1
c           write(6,*) itp
            temp(itp) = temper
            read (5,*) pres(itp),dkin(itp)
c           write (6,*) pres(itp),dkin(itp)
            pin(itp) = float(itp)
            sig(itp) = 1.0d0
            dkin(itp) = log10(dkin(itp))
         enddo
         if (it.eq.1) np1=np
      enddo
      ntp=itp
      dkltlp=10**dkin(1)
      dklthp=10**dkin(np1)
      dkhtlp=10**dkin(ntp-np+1)
      dkhthp=10**dkin(ntp)
      tmax = temp(ntp)
      tmin = temp(1)
      pminlt = pres(1)
      pminht = pres(ntp-np+1)
      dk0lt=dkltlp/pminlt
      dk0ht=dkhtlp/pminht
      if (iinp.eq.0) then
         a0(2)=0.0d0
         a0(5)=0.0d0
         a0(7)=0.7d0
         a0(8)=200.d0
         a0(9)=300.d0
         a0(10)=400.d0 
         a0(3) = (log(dk0ht)-log(dk0lt))/(1.0d0/tmin-1.0d0/tmax)
         a0(1) = log(dk0lt*exp(a0(3)/tmin))
         a0(6) = (log(dkhthp)-log(dklthp))/(1.0d0/tmin-1.0d0/tmax)
         a0(4) = log(dklthp*exp(a0(6)/tmin))
      endif
      write (6,102) exp(a0(1)),a0(2),a0(3)
      write (6,102) exp(a0(4)),a0(5),a0(6)
      write (6,102) exp(-a0(7)**2),(abs(a0(ii)),ii=8,10)
      if (ilist(9).ne.1) a0(9)=1.0d20
      nca=300
      alamda = -1.0
      iit = 0
      chisq = 0.0
      chisqo = 0.0
      iit = 0
  100 continue
      if (((chisq-chisqo).lt.0.0d0).and.
     $ (((chisqo-chisq).lt.converg1).or.
     $ ((chisqo-chisq)/chisq.lt.converg2))) go to 200
      chisqo=chisq
      call mrqmin(pin,dkin,sig,ntp,a0,ilist,ma,covar,alpha,
     $ nca,chisq,funcs,alamda)
      iit = iit + 1
      write (6,*) 'results for iteration',iit
      write (6,101) chisq,alamda
 101  format (1x,'chi ',g12.5,'lamda ',g12.5)
      write (6,*) 'params'
      write (6,102) exp(a0(1)),a0(2),a0(3)
      write (6,102) exp(a0(4)),a0(5),a0(6)
      write (6,102) exp(-a0(7)**2),(abs(a0(ii)),ii=8,10)
c     write (6,101) (a0(ii),ii=1,ma),chisq,alamda
 102  format (1x,4g12.5)
      do ip = 1 , ntp
         tempi = temp(ip)
         presi = pres(ip)
         do ii = 1 , ma
            if ((ii.eq.8).or.(ii.eq.9).or.(ii.eq.10)) then
               a(ii) = abs(a0(ii))
            else
               if (ii.eq.7) then
                  a(ii)=exp(-a0(ii)**2)
               else
                  a(ii)=a0(ii)
               endif
            endif
         enddo

c calculate reference fitted k

         dk0l = a(1)+a(2)*log(tempi)-a(3)/tempi
         dk0 = exp(dk0l)
         dkinfl = a(4)+a(5)*log(tempi)-a(6)/tempi
         dkinf = exp(dkinfl)
         fcent = (1.0d0-a(7))*exp(-tempi/a(10))+a(7)*exp(-tempi/a(8))+
     $    exp(-a(9)/tempi)
         pred = dk0*presi/dkinf
         dklind = dkinf*pred/(1.0d0+pred)
         cc = -0.4d0 - 0.67d0*dlog10(fcent)
         dn = 0.75d0 - 1.27d0*dlog10(fcent)
         tmp1 = dlog10(pred) + cc
         tmp2 = tmp1/(dn-0.14*tmp1)
         tmp3 = 1.0d0/(1.0d0+tmp2**2)
         ff = 10**(tmp3*dlog10(fcent))
         dkfit = dklind*ff
         dklow = dk0*presi
         write (6,111) tempi,presi,10**(dkin(ip)),dkfit,dklind,
     $    dklow,dkinf,fcent,(10**dkin(ip)-dkfit)*100.d0/dkfit
 111     format (1x,10g11.4)
      enddo
      go to 100
  200 continue

      stop
      end



      subroutine funcs(p,a0,dkfit,dkda,ma)

      implicit real*8(a-h,o-z)
      dimension a(ma),a0(ma),dkda(ma)
      dimension dkfitp(ma),dkfitm(ma),astep(ma)
      common /tp/ temp(300),pres(300)

      itp = idint(p+0.5d0)
      tempi = temp(itp)
      presi = pres(itp)
      do ii = 1 , ma
         if ((ii.eq.8).or.(ii.eq.9).or.(ii.eq.10)) then
            a(ii) = abs(a0(ii))
         else
c           if ((ii.eq.7).and.(abs(a0(ii)).gt.1.0d0)) then
            if (ii.eq.7) then
               a(ii)=exp(-a0(ii)**2)
            else
               a(ii)=a0(ii)
            endif
         endif
      enddo

c First calculate reference fitted k

      dk0l = a(1)+a(2)*log(tempi)-a(3)/tempi
      dk0 = exp(dk0l)
      dkinfl = a(4)+a(5)*log(tempi)-a(6)/tempi
      dkinf = exp(dkinfl)
      fcent = (1.0d0-a(7))*exp(-tempi/a(10))+a(7)*exp(-tempi/a(8))+
     $ exp(-a(9)/tempi)
      pred = dk0*presi/dkinf
      dklind = dkinf*pred/(1.0d0+pred)
      cc = -0.4d0 - 0.67d0*dlog10(fcent)
      dn = 0.75d0 - 1.27d0*dlog10(fcent)
      tmp1 = dlog10(pred) + cc
      tmp2 = tmp1/(dn-0.14*tmp1)
      tmp3 = 1.0d0/(1.0d0+tmp2**2)
      ff = 10**(tmp3*dlog10(fcent))
      dkfit = dlog10(dklind*ff)
c     write (6,*) 'dkfit test',dklind,ff

c Now calculate derivatives with respect to paramaters

      do ii = 1 , ma
         astep(ii) = 0.001d0*a0(ii)
         if ((ii.eq.2).or.(ii.eq.5)) astep(ii)=max(0.001d0,astep(ii))
         if ((ii.eq.8).or.(ii.eq.9).or.(ii.eq.10)) then
            a(ii) = abs(a0(ii)+astep(ii))
         else
c           if ((ii.eq.7).and.(abs(a0(ii)+astep(ii)).gt.1.0d0)) then
            if (ii.eq.7) then
               a(ii)=exp(-(a0(ii)+astep(ii))**2)
            else
               a(ii)=a0(ii)+astep(ii)
            endif
         endif
         dk0l = a(1)+a(2)*log(tempi)-a(3)/tempi
         dk0 = exp(dk0l)
         dkinfl = a(4)+a(5)*log(tempi)-a(6)/tempi
         dkinf = exp(dkinfl)
         fcent = (1.0d0-a(7))*exp(-tempi/a(10))+a(7)*exp(-tempi/a(8))+
     $    exp(-a(9)/tempi)
         pred = dk0*presi/dkinf
         dklind = dkinf*pred/(1.0d0+pred)
         cc = -0.4d0 - 0.67d0*dlog10(fcent)
         dn = 0.75d0 - 1.27d0*dlog10(fcent)
         tmp1 = dlog10(pred) + cc
         tmp2 = tmp1/(dn-0.14*tmp1)
         tmp3 = 1.0d0/(1.0d0+tmp2**2)
         ff = 10**(tmp3*dlog10(fcent))
         dkfitp(ii) = dlog10(dklind*ff)
         if ((ii.eq.8).or.(ii.eq.9).or.(ii.eq.10)) then
            a(ii) = abs(a0(ii))
         else
            if (ii.eq.7) then
               a(ii)=exp(-a0(ii)**2)
            else
               a(ii)=a0(ii)
            endif
         endif
c     write (6,*) 'dkfitp test',dklind,ff,dkfitp(ii),ii,astep(ii)
      enddo

      do ii = 1 , ma
         if ((ii.eq.8).or.(ii.eq.9).or.(ii.eq.10)) then
            a(ii) = abs(a0(ii)-astep(ii))
         else
            if (ii.eq.7) then
               a(ii)=exp(-(a0(ii)-astep(ii))**2)
            else
               a(ii)=a0(ii)-astep(ii)
            endif
         endif
         dk0l = a(1)+a(2)*log(tempi)-a(3)/tempi
         dk0 = exp(dk0l)
         dkinfl = a(4)+a(5)*log(tempi)-a(6)/tempi
         dkinf = exp(dkinfl)
         fcent = (1.0d0-a(7))*exp(-tempi/a(10))+a(7)*exp(-tempi/a(8))+
     $    exp(-a(9)/tempi)
         pred = dk0*presi/dkinf
         dklind = dkinf*pred/(1.0d0+pred)
         cc = -0.4d0 - 0.67d0*dlog10(fcent)
         dn = 0.75d0 - 1.27d0*dlog10(fcent)
         tmp1 = dlog10(pred) + cc
         tmp2 = tmp1/(dn-0.14*tmp1)
         tmp3 = 1.0d0/(1.0d0+tmp2**2)
         ff = 10**(tmp3*dlog10(fcent))
         dkfitm(ii) = dlog10(dklind*ff)
         if ((ii.eq.8).or.(ii.eq.9).or.(ii.eq.10)) then
            a(ii) = abs(a0(ii))
         else
            if (ii.eq.7) then
               a(ii)=exp(-a0(ii)**2)
            else
               a(ii)=a0(ii)
            endif
         endif
c     write (6,*) 'dkfitm test',dklind,ff,dkfitm(ii),ii,astep(ii)
      enddo
      do ii = 1 , ma
c     write (7,*) 'dkda test',ii,dkfitp(ii),dkfitm(ii),astep(ii)
         dkda(ii) = (dkfitp(ii)-dkfitm(ii))/(2.0d0*astep(ii))
      enddo
c     write (7,111) p,10**dkfit,(a(i),i=1,ma),(dkda(i),i=1,ma)
c111  format (1x,'func test',20g12.5)
c     write (7,112) (dkfitp(ii),ii=1,ma)
c     write (7,112) (dkfitm(ii),ii=1,ma)
 112  format (1x,6g12.5)
     
      return
      end 



      SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,
     *funcs,alamda)
      INTEGER ma,nca,ndata,ia(ma),MMAX
      REAL*8 alamda,chisq,funcs,a(ma),alpha(nca,nca),covar(nca,nca),
     *sig(ndata),x(ndata),y(ndata)
      PARAMETER (MMAX=20)
CU    USES covsrt,gaussj,mrqcof
      INTEGER j,k,l,mfit
      REAL*8 ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      SAVE ochisq,atry,beta,da,mfit
      if(alamda.lt.0.)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).ne.0) mfit=mfit+1
11      continue
        alamda=0.001
        call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,funcs)
        ochisq=chisq
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif
      do 14 j=1,mfit
        do 13 k=1,mfit
          covar(j,k)=alpha(j,k)
13      continue
        covar(j,j)=alpha(j,j)*(1.+alamda)
        da(j)=beta(j)
14    continue
      call gaussj(covar,mfit,nca,da,1,1)
      if(alamda.eq.0.)then
        call covsrt(covar,nca,ma,ia,mfit)
        return
      endif
      j=0
      do 15 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          atry(l)=a(l)+da(j)
        endif
15    continue
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,funcs)
      if(chisq.lt.ochisq)then
        alamda=0.1*alamda
        ochisq=chisq
        do 17 j=1,mfit
          do 16 k=1,mfit
            alpha(j,k)=covar(j,k)
16        continue
          beta(j)=da(j)
17      continue
        do 18 l=1,ma
          a(l)=atry(l)
18      continue
      else
        alamda=10.*alamda
        chisq=ochisq
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .


      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,
     *funcs)
      INTEGER ma,nalp,ndata,ia(ma),MMAX
      REAL*8 chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata),
     *y(ndata)
      EXTERNAL funcs
      PARAMETER (MMAX=20)
      INTEGER mfit,i,j,k,l,m
      REAL*8 dy,sig2i,wt,ymod,dyda(MMAX)
      mfit=0
      do 11 j=1,ma
        if (ia(j).ne.0) mfit=mfit+1
11    continue
      do 13 j=1,mfit
        do 12 k=1,j
          alpha(j,k)=0.
12      continue
        beta(j)=0.
13    continue
      chisq=0.
      do 16 i=1,ndata
        call funcs(x(i),a,ymod,dyda,ma)
        sig2i=1./(sig(i)*sig(i))
        dy=y(i)-ymod
        j=0
        do 15 l=1,ma
          if(ia(l).ne.0) then
            j=j+1
            wt=dyda(l)*sig2i
            k=0
            do 14 m=1,l
              if(ia(m).ne.0) then
                k=k+1
                alpha(j,k)=alpha(j,k)+wt*dyda(m)
              endif
14          continue
            beta(j)=beta(j)+dy*wt
          endif
15      continue
        chisq=chisq+dy*dy*sig2i
16    continue
      do 18 j=2,mfit
        do 17 k=1,j-1
          alpha(k,j)=alpha(j,k)
17      continue
18    continue
      return
      END


      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL*8 a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL*8 big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                write(6, *) 'singular matrix in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) write(6, *) 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .

      SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
      INTEGER ma,mfit,npc,ia(ma)
      REAL*8 covar(npc,npc)
      INTEGER i,j,k
      REAL*8 swap
      do 12 i=mfit+1,ma
        do 11 j=1,i
          covar(i,j)=0.
          covar(j,i)=0.
11      continue
12    continue
      k=mfit
      do 15 j=ma,1,-1
        if(ia(j).ne.0)then
          do 13 i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
13        continue
          do 14 i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
14        continue
          k=k-1
        endif
15    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .

