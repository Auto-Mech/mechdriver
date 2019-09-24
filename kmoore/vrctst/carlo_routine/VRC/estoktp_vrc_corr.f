
      double precision function estoktp_vrc_corr(natoms,x,rparm,iparm)
c
c     natoms = number of atoms (not used)
c     x = cartesian coordinates in au
c     e = energy correction in au
c     rparm (not used)
c
      implicit double precision (a-h,o-z)
      parameter (natommx=100)
      parameter (npointsmx=100)
      dimension x(3,natommx)
      dimension iparm(4)
      dimension rinp(npointsmx),dv(npointsmx),dv20(npointsmx)
      data dvp1,dvpn / 1.0d40,1.0d40 /
c      data rdistmin,rdistmax / 1.6d0, 10.0d0 /

      if(iparm(1).eq.1)then
         estoktp_vrc_corr = 0.
         goto 100
      endif
      open(unit=10,file='../../pot.dat',status='old',action='read')
      read(10,*)na,nb
      read(10,*)npoints
      do j=1,npointsmx
         rinp(j)=0.
         dv(j)=0.
         dv20(j)=0.
      enddo
      do j=1,npoints
         read(10,*)rinp(j),dv(j)
      enddo
      close(10)
      rdistmin=rinp(1)
      rdistmax=rinp(npoints)
      do j=1,npoints
         if(rinp(j).lt.rdistmin)rdistmin=rinp(j)
         if(rinp(j).gt.rdistmax)rdistmax=rinp(j)
      enddo
     
      rdist = dsqrt( (x(1,nb)-x(1,na))**2 + (x(2,nb)-x(2,na))**2 + 
     x             (x(3,nb)-x(3,na))**2 )

      rdist = rdist*0.529167
      delmlt = 1.0d0
      if(rdist.le.rdistmin) rdist = rdistmin
      if(rdist.ge.rdistmax) then
        delmlt = exp(-2.0d0*(rdist-rdistmax))
        rdist=rdistmax
      endif

      call spline(rinp,dv,npoints,dvp1,dvpn,dv20)
      call splint(rinp,dv,dv20,npoints,rdist,estoktp_vrc_corr)

      estoktp_vrc_corr = estoktp_vrc_corr*delmlt/627.5d0
      open(unit=10,file='test.dat',status='unknown')
      write (10,*) 'test',rdist,estoktp_vrc_corr*627.5
      close(10)

 100  continue

      return 
      end



