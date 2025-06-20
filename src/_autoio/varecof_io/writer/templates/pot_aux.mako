
      subroutine cross(a,b,c)
 
      implicit real*8(a-h,o-z)
      dimension a(3),b(3),c(3)
 
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
 
      return
      end
 

      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end


      subroutine dgtsl(n,c,d,e,b,info)
      integer n,info
      double precision c(1),d(1),e(1),b(1)
c
c     dgtsl given a general tridiagonal matrix and a right hand
c     side will find the solution.
c
c     on entry
c
c        n       integer
c                is the order of the tridiagonal matrix.
c
c        c       double precision(n)
c                is the subdiagonal of the tridiagonal matrix.
c                c(2) through c(n) should contain the subdiagonal.
c                on output c is destroyed.
c
c        d       double precision(n)
c                is the diagonal of the tridiagonal matrix.
c                on output d is destroyed.
c
c        e       double precision(n)
c                is the superdiagonal of the tridiagonal matrix.
c                e(1) through e(n-1) should contain the superdiagonal.
c                on output e is destroyed.
c
c        b       double precision(n)
c                is the right hand side vector.
c
c     on return
c
c        b       is the solution vector.
c
c        info    integer
c                = 0 normal value.
c                = k if the k-th element of the diagonal becomes
c                    exactly zero.  the subroutine returns when
c                    this is detected.
c
c     linpack. this version dated 08/14/78 .
c     jack dongarra, argonne national laboratory.
c
c     no externals
c     fortran dabs
c
c     internal variables
c
      integer k,kb,kp1,nm1,nm2
      double precision t
c     begin block permitting ...exits to 100
c
         info = 0
         c(1) = d(1)
         nm1 = n - 1
         if (nm1 .lt. 1) go to 40
            d(1) = e(1)
            e(1) = 0.0d0
            e(n) = 0.0d0
c
            do 30 k = 1, nm1
               kp1 = k + 1
c
c              find the largest of the two rows
c
               if (dabs(c(kp1)) .lt. dabs(c(k))) go to 10
c
c                 interchange row
c
                  t = c(kp1)
                  c(kp1) = c(k)
                  c(k) = t
                  t = d(kp1)
                  d(kp1) = d(k)
                  d(k) = t
                  t = e(kp1)
                  e(kp1) = e(k)
                  e(k) = t
                  t = b(kp1)
                  b(kp1) = b(k)
                  b(k) = t
   10          continue
c
c              zero elements
c
               if (c(k) .ne. 0.0d0) go to 20
                  info = k
c     ............exit
                  go to 100
   20          continue
               t = -c(kp1)/c(k)
               c(kp1) = d(kp1) + t*d(k)
               d(kp1) = e(kp1) + t*e(k)
               e(kp1) = 0.0d0
               b(kp1) = b(kp1) + t*b(k)
   30       continue
   40    continue
         if (c(n) .ne. 0.0d0) go to 50
            info = n
         go to 90
   50    continue
c
c           back solve
c
            nm2 = n - 2
            b(n) = b(n)/c(n)
            if (n .eq. 1) go to 80
               b(nm1) = (b(nm1) - d(nm1)*b(n))/c(nm1)
               if (nm2 .lt. 1) go to 70
               do 60 kb = 1, nm2
                  k = nm2 - kb + 1
                  b(k) = (b(k) - d(k)*b(k+1) - e(k)*b(k+2))/c(k)
   60          continue
   70          continue
   80       continue
   90    continue
  100 continue
c
      return
      end

      subroutine hunt(xx,n,x,jlo)
      implicit real*8 (a-h,o-z)
      dimension xx(n)
      logical ascnd

      ascnd = xx(n) .gt. xx(1)
      if( jlo .le. 0  .or. jlo .gt. n) then
         jlo = 0
         jhi=n+1
         go to 3
      endif
      inc=1

      if( x.ge. xx(jlo) .eqv. ascnd) then
 1       jhi = jlo+inc
         if(jhi .gt. n) then
           jhi=n+1
         else if( x .ge. xx(jhi) .eqv. ascnd) then
           jlo=jhi
           inc=inc+inc
           go to 1
         endif
      else
         jhi=jlo
 2       jlo=jhi-inc
         if(jlo .lt.1) then
            jlo=0
         else if( x .lt. xx(jlo) .eqv. ascnd) then
            jhi=jlo
            inc=inc+inc
            go to 2
         endif
      endif

 3    if(jhi -jlo .eq.1) return
      jm=(jhi+jlo)/2
      if( x .gt. xx(jm) .eqv. ascnd) then
         jlo = jm
      else
         jhi = jm
      endif
      go to 3
      end

      subroutine seval(c,px,py,pz,fi)
      implicit real*8 (a-h,o-z)
      dimension c(4,4,4),px(4),py(4),pz(4)
      dimension temp1(4,4),temp2(4)
c
      do 50 i=1,4
      do 50 j=1,4
   50 temp1(i,j)=ddot(4,c(1,i,j),1,pz,1)
      do 60 i=1,4
   60 temp2(i)=ddot(4,temp1(1,i),1,py,1)
      fi=ddot(4,temp2,1,px,1)
c
      return
      end


      subroutine spline(x,y,n,yp1,ypn,y2)
      integer n,nmax
      real*8 yp1,ypn,x(n),y(n),y2(n)
      parameter (nmax=500)
      integer i,k
      real*8 p,qn,sig,un,u(nmax)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end
c  (c) copr. 1986-92 numerical recipes software .


      subroutine splint(xa,ya,y2a,n,x,y)
      integer n
      real*8 x,y,xa(n),y2a(n),ya(n)
      integer k,khi,klo
      real*8 a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
c      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      end
c  (c) copr. 1986-92 numerical recipes software .


      double precision function plgndr(l,m,x)
      implicit real*8 (a-h,o-z)
c
c     from numerical recipes,
c
c     computes associated legendre polynomials by recurrence
c
      if(m.lt.0 .or. m.gt.l .or. abs(x) .gt. 1.0) then
        stop 'error in plgndr'
      endif
      pmm=1.0d0
      if(m .gt. 0) then
        somx2=sqrt((1.0d0-x)*(1.0d0+x))
        fact=1.0d0
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.0d0
 11     continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
 12       continue
          plgndr=pll
        endif
      endif

      return
      end

