      real*8 function dummy_corr(na,x,rparm,iparm)
c
c     na = number of atoms
c     x = cartesian coordinates
c     e = energy correction in au
c
      implicit real*8 (a-h,o-z)
      dimension x(3,1)
      dummy_corr = 0.0
      return
      end
