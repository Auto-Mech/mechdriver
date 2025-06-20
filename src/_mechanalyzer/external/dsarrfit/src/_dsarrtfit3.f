      Program dsarrfit

c Program to fit both
c (i) A,n,E in k=A T^n exp(-E/T)
c (ii) A1,n1,E1,A2,n2,E2 in k=sum_i Ai T^ni exp(-Ei/T)
c with the parameters specified in that order.

      implicit real*8 (a-h,o-z)
      parameter (mas=3,ma=6)
      dimension temp(100),tempp(100),dkin(100),sig(100)
      dimension a0(ma),a(ma),lista(ma),dkda(ma)
      dimension ilist(ma)
      dimension covar(100,100),alpha(100,100)
      character*2 ftype
      external funcs,funcd

      OPEN(UNIT=5,STATUS='OLD',FILE='arrfit.dat')

      read (5,*) converg1,converg2
      read (5,*) ftype
      read (5,*) (ilist(ii),ii=1,ma)
      read (5,*) iinp
      if (iinp.eq.1) then
         read (5,*) (a0(ii),ii=1,mas)
         a0(1) = log(a0(1))-a0(2)*log(298.)
      elseif (iinp.eq.2) then
         read (5,*) (a0(ii),ii=1,mas)
         read (5,*) (a0(ii+3),ii=1,mas)
         a0(1) = log(a0(1))-a0(2)*log(298.)
         a0(4) = log(a0(4))-a0(5)*log(298.)
      endif
      tmin = 1.0d30
      tmax = 0.0d0
      read (5,*) nt,nttot
      do it = 1 , nttot
         read (5,*) temp(it),dkin(it)
         sig(it) = 1.0d0
         dkin(it) = log10(dkin(it))
         if (it.gt.nt) go to 50
         if (temp(it).lt.tmin) then
            tmin=temp(it)
            dklt=10**dkin(it)
         endif
         if (temp(it).gt.tmax) then
            tmax=temp(it)
            dkht=10**dkin(it)
         endif
  50     continue
      enddo

c start with single modified arrhenius  **************************************

      if (ftype.eq.'s' .or. ftype.eq.'sd') then ! start single fit
      OPEN(UNIT=6,STATUS='unknown',FILE='sgl_arrfit.out')
c     generate initial guess if one wasn't given
      if (iinp.eq.0) then
         a0(3)=(log(dkht)-log(dklt))/(1.0d0/tmin-1.0d0/tmax)
         a0(1)=log(dklt*exp(a0(3)/tmin))
         a0(2)=0.0d0
      endif
c     print out initial guess comparison
      write (6,*) 'initial guess'
      write (6,102) exp(a0(1))*(298.**a0(2)),(a0(ii),ii=2,mas)
      do it = 1 , nttot
         tempi = temp(it)
         dkfitl = a0(1)+a0(2)*log(tempi)-a0(3)/tempi
         dkfit = exp(dkfitl)
         dki=10**(dkin(it))
         err = abs(dkfit-dki)*100./dki
         write (6,111) tempi,10**(dkin(it)),dkfit,err
      enddo
      nca=100
      mfit = 0
      do ii = 1 ,mas 
         if (ilist(ii).eq.1) then
            mfit = mfit + 1
            lista(mfit)=ii
         endif
      enddo 
      alamda = -1.0
      iit = 0
      chisq = 0.0
      chisqo = 0.0
      iit = 0
      iconverge = 0
  100 continue
      if (((chisq-chisqo).lt.0.0d0).and.
     $ (((chisqo-chisq).lt.converg1).or.
     $ ((chisqo-chisq)/chisq.lt.converg2))) go to 200
	iconverge = iconverge + 1
	if (iconverge .gt. 1000) go to 200
      chisqo=chisq
      call mrqmin(temp,dkin,sig,nt,a0,ilist,mas,covar,alpha,
     $ nca,chisq,funcs,alamda)
      iit = iit + 1
      write (6,*) 'results for iteration',iit
      write (6,101) chisq,alamda
 101  format (1x,'chi ',g12.5,' lamda ',g12.5)
      write (6,*) 'params'
c     write (6,102) exp(a0(1))*(298.**a0(2)),(a0(ii),ii=2,mas)
      write (6,102) exp(a0(1)),(a0(ii),ii=2,mas)
 102  format (1x,5(2x,ES13.5E3))
c 102  format (1x,5(2x,g13.5))
      write (6,102) 6.0221e23*exp(a0(1)),a0(2),a0(3)*1.987
      do it = 1 , nttot
         tempi = temp(it)
c calculate reference fitted k
         dkfitl = a0(1)+a0(2)*log(tempi)-a0(3)/tempi
         dkfit = exp(dkfitl)
         dki=10**(dkin(it))
         err = abs(dkfit-dki)*100./dki
         write (6,111) tempi,10**(dkin(it)),dkfit,err
 111     format (1x,5g12.5)
      enddo
      go to 100
  200 continue
c     Write the converged results
      write (6,*)
      write (6,*)
      write (6,*) "Converged Params for single modified Arrhenius"
      write (6,101) chisq,alamda
      write (6,*) 'params'
c      write (6,1002) a0(1)*(298.**a0(2)),a0(2),a0(3),
c     $ a0(4)*(298.**a0(5)),a0(5),a0(6)
c      write (6,1002) dlog10(a0(1)),a0(2),a0(3)*1.987,
c     $ dlog10(a0(4)),a0(5),a0(6)*1.987
      write (6,*) '   A              n              Ea/R'
      write (6,102) exp(a0(1)),(a0(ii),ii=2,mas)
      write (6,*) '   A*Navo         n              Ea (kcal/mol)'
      write (6,102) 6.0221e23*exp(a0(1)),a0(2),a0(3)*1.987
      do it = 1 , nttot
         tempi = temp(it)
c calculate reference fitted k
         dkfitl1 = a0(1)+a0(2)*log(tempi)-a0(3)/tempi
         dkfit = exp(dkfitl1)
         dki=10**(dkin(it))
         err = abs(dkfit-dki)*100./dki
         write (6,111) tempi,10**(dkin(it)),dkfit,err
      enddo
      if (nttot.lt.4) then
         go to 8000
      endif
      endif ! end single fit
c
c now repeat for sum of two modified arrhenius    *****************************
c
      if (ftype.eq.'d' .or. ftype.eq.'sd') then ! start double fit
      OPEN(UNIT=7,STATUS='unknown',FILE='dbl_arrfit.out')
      if (iinp.eq.0 .or. iinp.eq.1 .or. ftype.eq.'sd') then
          a1inp = exp(a0(1))*(298.**a0(2))/2.0d0
          a0(2) = a0(2)-2.0d0
          a0(5) = a0(2)+4.0d0
          a0(6) = a0(3)
          a0(1) = a1inp/(298.**a0(2))
          a0(4) = a1inp/(298.**a0(5))
      endif
      write (7,*) 'initial guesses'
      write (7,1002) a0(1)*(298.**a0(2)),a0(2),a0(3),
     $ a0(4)*(298.**a0(5)),a0(5),a0(6)
      do it = 1 , nttot
         tempi = temp(it)
c calculate reference fitted k
         dkfitl1 = a0(2)*log(tempi)-a0(3)/tempi
         dkfitl2 = a0(5)*log(tempi)-a0(6)/tempi
         dkfit = a0(1)*exp(dkfitl1)+a0(4)*exp(dkfitl2)
         dki=10**(dkin(it))
         err = abs(dkfit-dki)*100./dki
         write (7,111) tempi,10**(dkin(it)),dkfit,err
      enddo
      nca=100
      mfit = 0
      do ii = 1 ,ma 
         if (ilist(ii).eq.1) then
            mfit = mfit + 1
            lista(mfit)=ii
         endif
      enddo 
      alamda = -1.0
      iit = 0
      chisq = 0.0
      chisqo = 0.0
      iit = 0
      iconverge = 0
 1000 continue
      if (((chisq-chisqo).lt.0.0d0).and.
     $ (((chisqo-chisq).lt.converg1).or.
     $ ((chisqo-chisq)/chisq.lt.converg2))) go to 2000
	iconverge = iconverge + 1
	if (iconverge .gt. 100000) go to 2000
      chisqo=chisq
      call mrqmin(temp,dkin,sig,nt,a0,ilist,ma,covar,alpha,
     $ nca,chisq,funcd,alamda)
      iit = iit + 1
      write (7,*) 'results for iteration',iit
      write (7,1001) chisq,alamda
 1001 format (1x,'chi ',g12.5,' lamda ',g12.5)
      write (7,*) 'params'
      write (7,1002) a0(1)*(298.**a0(2)),a0(2),a0(3),
     $ a0(4)*(298.**a0(5)),a0(5),a0(6)
 1002 format (1x,6(2x,ES16.5E3))
c 1002 format (1x,6(2x,g12.5))
      do it = 1 , nttot
         tempi = temp(it)
c calculate reference fitted k
         dkfitl1 = a0(2)*log(tempi)-a0(3)/tempi
         dkfitl2 = a0(5)*log(tempi)-a0(6)/tempi
         dkfit = a0(1)*exp(dkfitl1)+a0(4)*exp(dkfitl2)
         dki=10**(dkin(it))
         err = abs(dkfit-dki)*100./dki
         write (7,111) tempi,10**(dkin(it)),dkfit,err
      enddo
      go to 1000
 2000 continue
c     Write the converged results
      write (7,*)
      write (7,*)
      write (7,*) "Converged Params from double modified Arrhenius"
      write (7,1001) chisq,alamda
      write (7,*) 'params'
c      write (7,1002) a0(1)*(298.**a0(2)),a0(2),a0(3),
c     $ a0(4)*(298.**a0(5)),a0(5),a0(6)
c      write (7,1002) dlog10(a0(1)),a0(2),a0(3)*1.987,
c     $ dlog10(a0(4)),a0(5),a0(6)*1.987
      write (7,*) 'A              n              Ea/R'
      write (7,1002) a0(1),a0(2),a0(3),
     $ a0(4),a0(5),a0(6)
      write (7,*) 'A*Navo         n              Ea (kcal/mol)'
      write (7,1002) 6.0221e23*a(1),a0(2),a0(3)*1.987,
     $ 6.0221e23*a0(4),a0(5),a0(6)*1.987
      do it = 1 , nttot
         tempi = temp(it)
c calculate reference fitted k
         dkfitl1 = a0(2)*log(tempi)-a0(3)/tempi
         dkfitl2 = a0(5)*log(tempi)-a0(6)/tempi
         dkfit = a0(1)*exp(dkfitl1)+a0(4)*exp(dkfitl2)
         dki=10**(dkin(it))
         err = abs(dkfit-dki)*100./dki
         write (7,111) tempi,10**(dkin(it)),dkfit,err
      enddo
      endif  ! end of double fit
 8000 continue
      stop
      end
cccccccccc
cccccccccc
      subroutine funcs(temp,a0,dkfit,dkda,ma)
      implicit real*8(a-h,o-z)
      dimension a(ma),a0(ma),dkda(ma)
      dimension dkfitp(ma),dkfitm(ma),astep(ma)
c First calculate reference fitted k
      do ii = 1 , ma
         a(ii)=a0(ii)
      enddo
      dkfitl = a0(1)+a0(2)*log(temp)-a0(3)/temp
      dkfit = log10(exp(dkfitl))
c Now calculate derivatives with respect to paramaters
      do ii = 1 , ma
         astep(ii) = 0.001d0*a0(ii)
         if (ii.eq.2) astep(ii)=max(0.001d0,astep(ii))
         a(ii)=a0(ii)+astep(ii)
         dkfitl = a(1)+a(2)*log(temp)-a(3)/temp
         dkfitp(ii) = log10(exp(dkfitl))
         a(ii)=a0(ii)
      enddo
      do ii = 1 , ma
         if (ii.eq.2) astep(ii)=max(0.001d0,astep(ii))
         a(ii)=a0(ii)-astep(ii)
         dkfitl = a(1)+a(2)*log(temp)-a(3)/temp
         dkfitm(ii) = log10(exp(dkfitl))
         a(ii)=a0(ii)
      enddo
      do ii = 1 , ma
c     write (7,*) 'dkda test',ii,dkfitp(ii),dkfitm(ii),astep(ii)
         dkda(ii) = (dkfitp(ii)-dkfitm(ii))/(2.0d0*astep(ii))
      enddo
c     write (7,111) p,10**dkfit,(a(i),i=1,ma),(dkda(i),i=1,ma)
 111  format (1x,'func test',20g12.5)
c     write (7,112) (dkfitp(ii),ii=1,ma)
c     write (7,112) (dkfitm(ii),ii=1,ma)
 112  format (1x,6g12.5)
     
      return
      end 
cccccccccc
cccccccccc
      subroutine funcd(temp,a0,dkfit,dkda,ma)
      implicit real*8(a-h,o-z)
      dimension a(ma),a0(ma),dkda(ma)
      dimension dkfitp(ma),dkfitm(ma),astep(ma)
c First calculate reference fitted k
      do ii = 1 , ma
         a(ii)=a0(ii)
      enddo
      dkfitl1 = a0(2)*log(temp)-a0(3)/temp
      dkfitl2 = a0(5)*log(temp)-a0(6)/temp
      if ((dkfitl1.gt.300.).or.(dkfitl2.gt.300.)) then
         dkfit = 1.0d100
      else
         dktmp = a0(1)*exp(dkfitl1)+a0(4)*exp(dkfitl2)
         dkfit = dktmp*log10(abs(dktmp))/abs(dktmp)
C        dkfit = log10(abs(a0(1)*exp(dkfitl1)+a0(4)*exp(dkfitl2)))
      endif
c Now calculate derivatives with respect to paramaters
      do ii = 1 , ma
         astep(ii) = 0.001d0*a0(ii)
         if ((ii.eq.2).or.(ii.eq.5)) astep(ii)=max(0.001d0,astep(ii))
         a(ii)=a0(ii)+astep(ii)
         dkfitl1 = a(2)*log(temp)-a(3)/temp
         dkfitl2 = a(5)*log(temp)-a(6)/temp
         dkfitp(ii) = log10(abs(a(1)*exp(dkfitl1)+a(4)*exp(dkfitl2)))
         a(ii)=a0(ii)
      enddo
      do ii = 1 , ma
         if ((ii.eq.2).or.(ii.eq.5)) astep(ii)=max(0.001d0,astep(ii))
         a(ii)=a0(ii)-astep(ii)
         dkfitl1 = a(2)*log(temp)-a(3)/temp
         dkfitl2 = a(5)*log(temp)-a(6)/temp
         dkfitm(ii) = log10(abs(a(1)*exp(dkfitl1)+a(4)*exp(dkfitl2)))
         a(ii)=a0(ii)
      enddo
      do ii = 1 , ma
c     write (7,*) 'dkda test',ii,dkfitp(ii),dkfitm(ii),astep(ii)
         dkda(ii) = (dkfitp(ii)-dkfitm(ii))/(2.0d0*astep(ii))
      enddo
c     write (7,111) 10**dkfit,(a(i),i=1,ma),(dkda(i),i=1,ma)
 111  format (1x,'func test',20g12.5)
c     write (7,112) (dkfitp(ii),ii=1,ma)
c     write (7,112) (dkfitm(ii),ii=1,ma)
 112  format (1x,6g12.5)
     
      return
      end 
cccccccccc
cccccccccc
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
cccccccccc
cccccccccc
C  (C) Copr. 1986-92 Numerical Recipes Software .
      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,
     *funcs)
      INTEGER ma,nalp,ndata,ia(ma),MMAX
      REAL*8 chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata),
     *y(ndata)
c     EXTERNAL funcs
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
		 stop
 !               pause 'singular matrix in gaussj'
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
        if (a(icol,icol).eq.0.) stop !pause 'singular matrix in gaussj'
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
