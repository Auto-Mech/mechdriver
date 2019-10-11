      subroutine bmatrix(ispecies,ifile)

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

      dimension tau(ntaumx),taumn(ntaumx),taumx(ntaumx),freq(nmdmx)
      dimension coord(natommx,ndim),xint(3*natommx),xinti(3*natommx),
     $ abcrot(ndim),xintt(3*natommx)

      dimension xintp(3*natommx),xintn(3*natommx)
      dimension xintpp(3*natommx),xintpn(3*natommx)
      dimension xintnp(3*natommx),xintnn(3*natommx)

      dimension bmat(3*natommx,3*natommx)
      dimension cmat(3*natommx,3*natommx,3*natommx)

      dimension coox(natommx),cooy(natommx),cooz(natommx)
      dimension cooxn(natommx),cooyn(natommx),coozn(natommx)
      dimension cooxp(natommx),cooyp(natommx),coozp(natommx)

      dimension cooxpp(natommx),cooypp(natommx),coozpp(natommx)
      dimension cooxnp(natommx),cooynp(natommx),cooznp(natommx)
      dimension cooxpn(natommx),cooypn(natommx),coozpn(natommx)
      dimension cooxnn(natommx),cooynn(natommx),cooznn(natommx)

      dimension cooxyz(3*natommx)
      dimension ibconn(natommx),iaconn(natommx),idconn(natommx)
      dimension idummy(natommx),iangind(3*natommx)
      dimension ianginda(3*natommx),iangindd(3*natommx)
      dimension idihed(natommx)
      dimension bval(natommx),aval(natommx),dval(natommx)


      LOGICAL leof,lsec,ltit
      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      character*20 bname(natommx),anname(natommx),dname(natommx)
     $ ,atname(natommx),bconnt(natommx),aconnt(natommx),dconnt(natommx)

      character*20 filename
      character*20 stoichname
      character*60 command1

      character*60 atomlabel(natommx)
      character*2 aname
      character*1 aname1(natommx)
      character*30 cjunk
      character*30 angsub1,angsub2
      character*30 nameout,nameoutc
      character*30 gmem
      character*30 intcoor(3*natommx)
      character*30 intcoort(3*natommx)
      character*20 bislab(ntaumx)

      include 'filcomm.f'


c read natom, charge, spin, ntau, zmat, islinear

cc initialize vectors
      do j=1,ncoord 
         bname(j)=''
         anname(j)=''
         dname(j)=''
         bval(j)=0.
         aval(j)=0.
         dval(j)=0.
      enddo
      do j=1,natomt
         ibconn(j)=0
         iaconn(j)=0
         idconn(j)=0
      enddo
cc
      call read_zmat(atomlabel,natom,natomt,intcoor,bislab,ibconn,
     $ iaconn,idconn,bname,anname,dname,atname,idummy,isited,jsited,
     $ ksited,bconnt,aconnt,dconnt)
 
cc this does not account for dummy atoms defined with  given quantities

      do k=1,natomt
         do j=1,ncoord
            if(intcoor(j).eq.bname(k))bval(k)=xint(j)
            if(intcoor(j).eq.anname(k))aval(k)=xint(j)
            if(intcoor(j).eq.dname(k))dval(k)=xint(j)
         enddo
      enddo

cc all that follows is not necessary. we read the intial geometry
cc I leave (coomented) it as it can be useful sooner or later

cc   define starting xyz geometry. 
cc   convention: atom 1 is is 0 0 0
cc   atom 2 bd 0 0
cc   atom 3 on xy plane 

cc read the zmat from the level 1 output
      read(18,*)natomcar
      read(18,*)cjunk
      do j=1,natomcar
         read(18,*)cjunk,coox(j),cooy(j),cooz(j)
c         write(*,*)cjunk,coox(j),cooy(j),cooz(j)
      enddo
      close(18)

cc this is a check that we can compute the xyz matrix and convert back to internal
      ntau=0
      ifilu=1
      call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $    ,idconn,bname,anname,dname,atname,coox,cooy,cooz,xinti,tauopt,
     $     ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)
      if(xint(1).eq.0)then
         write(*,*)'failed in updating geometry'
         write(*,*)'probably error of g09, check output'
         write(*,*)'stopping now'
         stop
      endif

      open(unit=10,file='temp1.xyz',status='unknown')
      do j=1,ncoord
         write(10,*)intcoor(j),xint(j),xinti(j),abs(xint(j)-xinti(j))
      enddo
      close(10)

      do j=1,ncoord
         xint(j)=xinti(j)
      enddo

c      stop

      open(unit=10,file='temp2.xyz',status='unknown')
      write(10,*)natom
      write(10,*)'test'
      inda=0
      do j=1,natom
         aname1(j)=atname(j)
         if(idummy(j).ne.1)then
            inda=inda+1
      endif
      write(10,*)aname1(inda),coox(j),cooy(j),cooz(j)
      enddo
      close(10)

c      stop

cc now we update the coordinates using those etransformed in our procedure for consistentcy
cc with dimension of dihedral angles

      do j=1,ncoord
         xint(j)=xinti(j)
      enddo

cc here I inizialize indexes to check if int coor is an angle or a dihedral
cc

      do k=1,ncoord
         iangind(k)=0
         ianginda(k)=0
         iangindd(k)=0
      enddo

      do k=1,ncoord
         do j=1,ncoord
            if(intcoor(k).eq.anname(j))iangind(k)=1
            if(intcoor(k).eq.anname(j))ianginda(k)=1
            if(intcoor(k).eq.dname(j))iangind(k)=1
            if(intcoor(k).eq.dname(j))iangindd(k)=1
         enddo
      enddo

cc now we can compute the B matrix using central difference
cc Bik=dqi/dxk

      deltax=0.01
      deltay=0.01
      deltaz=0.01
      kind=0
      kaind=0
      ifilu=0

      do ij=1,ncoord
         do ik=1,3*natom
            bmat(ij,ik)=0.
         enddo
      enddo

      do k=1,natomt
         if(idummy(k).ne.1)then
            kind=kind+1
            kaind=kaind+1
            write(*,*)'k atom is ', k
            write(*,*)'kaind is ', kaind
            write(*,*)'kind is ', kind
            do ij=1,ncoord
               xintp(ij)=0.
            enddo

            do ik=1,natom
               cooxp(ik)=coox(ik)
               cooyp(ik)=cooy(ik)
               coozp(ik)=cooz(ik)
            enddo
            cooxp(kaind)=cooxp(kaind)+deltax

            call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $    ,idconn,bname,anname,dname,atname,cooxp,cooyp,coozp,xintp,
     $  tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)

            do ik=1,natom
               cooxn(ik)=coox(ik)
               cooyn(ik)=cooy(ik)
               coozn(ik)=cooz(ik)
            enddo
            cooxn(kaind)=cooxn(kaind)-deltax

            call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $    ,idconn,bname,anname,dname,atname,cooxn,cooyn,coozn,xintn,
     $  tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)

cc we can now compute the first component of the Bmatrix

            do i=1,ncoord
               if(abs(xintp(i)-xintn(i)).gt.300)then
                  if(xintn(i).lt.0)then
                     xintn(i)=xintn(i)+360.
                  else if(xintn(i).gt.0)then
                     xintn(i)=xintn(i)-360.
                  endif
                  if(abs(xintp(i)-xintn(i)).gt.300)then
                     write(*,*)'something did not work here'
                     write(*,*)'kind= ',kind
                     write(*,*)'intocoor= ',i
                     write(*,*)'xintp= ',xintp(i)
                     write(*,*)'xintn= ',xintn(i)
                     stop
                  endif
               endif
               bmat(i,kind)=(xintp(i)-xintn(i))/2./deltax
            enddo

cc we now replicate for y

cc first perturb y + deltay
            do ik=1,natom
               cooxp(ik)=coox(ik)
               cooyp(ik)=cooy(ik)
               coozp(ik)=cooz(ik)
            enddo
            cooyp(kaind)=cooyp(kaind)+deltay

            call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $    ,idconn,bname,anname,dname,atname,cooxp,cooyp,coozp,xintp,
     $  tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)

cc then perturb y - deltay
            do ik=1,natom
               cooxn(ik)=coox(ik)
               cooyn(ik)=cooy(ik)
               coozn(ik)=cooz(ik)
            enddo
            cooyn(kaind)=cooyn(kaind)-deltay

            call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $    ,idconn,bname,anname,dname,atname,cooxn,cooyn,coozn,xintn,
     $  tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)

cc we can now compute the second(y) component of the Bmatrix for atom k

            kind=kind+1
           do i=1,ncoord
               if(abs(xintp(i)-xintn(i)).gt.300)then
                  if(xintn(i).lt.0)then
                     xintn(i)=xintn(i)+360.
                  else if(xintn(i).gt.0)then
                     xintn(i)=xintn(i)-360.
                  endif
                  if(abs(xintp(i)-xintn(i)).gt.300)then
                     write(*,*)'something did not work here'
                     write(*,*)'kind= ',kind
                     write(*,*)'intocoor= ',i
                     write(*,*)'xintp= ',xintp(i)
                     write(*,*)'xintn= ',xintn(i)
                     stop
                  endif
               endif
               bmat(i,kind)=(xintp(i)-xintn(i))/2./deltay
            enddo

cc we now replicate for z

cc first perturb z + deltaz
            do ik=1,natom
               cooxp(ik)=coox(ik)
               cooyp(ik)=cooy(ik)
               coozp(ik)=cooz(ik)
            enddo
            coozp(kaind)=coozp(kaind)+deltaz

            call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $    ,idconn,bname,anname,dname,atname,cooxp,cooyp,coozp,xintp,
     $  tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)

cc then perturb z - deltaz
            do ik=1,natom
               cooxn(ik)=coox(ik)
               cooyn(ik)=cooy(ik)
               coozn(ik)=cooz(ik)
            enddo
            coozn(kaind)=coozn(kaind)-deltaz

            call update_zmat(natom,natomt,intcoor,bislab,ibconn,iaconn
     $    ,idconn,bname,anname,dname,atname,cooxn,cooyn,coozn,xintn,
     $  tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,dconnt,atomlabel,ifilu)

cc we can now compute the third(z) component of the Bmatrix for atom k

            kind=kind+1
            do i=1,ncoord
               if(abs(xintp(i)-xintn(i)).gt.300)then
                  if(xintn(i).lt.0)then
                     xintn(i)=xintn(i)+360.
                  else if(xintn(i).gt.0)then
                     xintn(i)=xintn(i)-360.
                  endif
                  if(abs(xintp(i)-xintn(i)).gt.300)then
                     write(*,*)'something did not work here'
                     write(*,*)'kind= ',kind
                     write(*,*)'intocoor= ',i
                     write(*,*)'xintp= ',xintp(i)
                     write(*,*)'xintn= ',xintn(i)
                     stop
                  endif
               endif
               bmat(i,kind)=(xintp(i)-xintn(i))/2./deltaz
            enddo
         endif
      enddo

      do k=1,ncoord
         do j=1,3*natom
            if(iangind(k).eq.1)bmat(k,j)=bmat(k,j)/180.*pigr
         enddo
      enddo

      if(nintcoord.ne.0)then

         do j=1,ncoord
            if(intcoor(j).eq.angsub1)ibmatsub1=j
            if(intcoor(j).eq.angsub2)ibmatsub2=j
         enddo

         iatsub1=0
         iatsub2=0
         do j=1,natomt
            if(anname(j).eq.angsub1)iatsub1=j
            if(dname(j).eq.angsub2)iatsub2=j
         enddo
         
c     now update iangsub1 bmat component
         iatom=1
         iinda=0
         bcentx=0.
         do j=1,3*natom
            iinda=iinda+1
            if(iinda.eq.4)then
               iinda=1
               iatom=iatom+1
            endif
            if((iatom.ne.iatsub1).and.iatom.ne.ibconn(iatsub1).and.
     +          iatom.ne.iaconn(iatsub1))then
               bmat(ibmatsub1,j)=0.
            endif
         enddo
         bmat(ibmatsub1,(ibconn(iatsub1)-1)*3+1)=
     +       -(bmat(ibmatsub1,(iatsub1-1)*3+1)+
     +         bmat(ibmatsub1,(iaconn(iatsub1)-1)*3+1))/2.
         bmat(ibmatsub1,(ibconn(iatsub1)-1)*3+2)=
     +       -(bmat(ibmatsub1,(iatsub1-1)*3+2)+
     +         bmat(ibmatsub1,(iaconn(iatsub1)-1)*3+2))/2.
         bmat(ibmatsub1,(ibconn(iatsub1)-1)*3+3)=
     +       -(bmat(ibmatsub1,(iatsub1-1)*3+3)+
     +         bmat(ibmatsub1,(iaconn(iatsub1)-1)*3+3))/2.

cc now determine displacement in orthogonal plane

cc  first get components of vector on plane
         vpx=coox(iatsub1)-coox(ibconn(iatsub1))
         vpy=cooy(iatsub1)-cooy(ibconn(iatsub1))
         vpz=cooz(iatsub1)-cooz(ibconn(iatsub1))
         vpnorm=sqrt(vpx*vpx+vpy*vpy+vpz*vpz)
         vpx=vpx/vpnorm
         vpy=vpy/vpnorm
         vpz=vpz/vpnorm
cc now update components of of 2nd  bmat linear angle
cc starting values
         bm1x=bmat(ibmatsub1,(iatsub1-1)*3+1)
         bm1y=bmat(ibmatsub1,(iatsub1-1)*3+2)
         bm1z=bmat(ibmatsub1,(iatsub1-1)*3+3)
         bm2x=bmat(ibmatsub1,(ibconn(iatsub1)-1)*3+1)
         bm2y=bmat(ibmatsub1,(ibconn(iatsub1)-1)*3+2)
         bm2z=bmat(ibmatsub1,(ibconn(iatsub1)-1)*3+3)
         bm3x= bmat(ibmatsub1,(iaconn(iatsub1)-1)*3+1)
         bm3y= bmat(ibmatsub1,(iaconn(iatsub1)-1)*3+2)
         bm3z= bmat(ibmatsub1,(iaconn(iatsub1)-1)*3+3)
cc project to perpendicular plane
         bmn1x=bm1y*vpz-bm1z*vpy
         bmn1y=-bm1x*vpz+bm1z*vpx
         bmn1z=bm1x*vpy-bm1y*vpx
         bmn2x=bm2y*vpz-bm2z*vpy
         bmn2y=-bm2x*vpz+bm2z*vpx
         bmn2z=bm2x*vpy-bm2y*vpx
         bmn3x=bm3y*vpz-bm3z*vpy
         bmn3y=-bm3x*vpz+bm3z*vpx
         bmn3z=bm3x*vpy-bm3y*vpx

cc update vector.
         do j=1,3*natom
            bmat(ibmatsub2,j)=0.
         enddo
         bmat(ibmatsub2,(iatsub1-1)*3+1)=bmn1x
         bmat(ibmatsub2,(iatsub1-1)*3+2)=bmn1y
         bmat(ibmatsub2,(iatsub1-1)*3+3)=bmn1z

         bmat(ibmatsub2,(ibconn(iatsub1)-1)*3+1)=bmn2x
         bmat(ibmatsub2,(ibconn(iatsub1)-1)*3+2)=bmn2y
         bmat(ibmatsub2,(ibconn(iatsub1)-1)*3+3)=bmn2z

         bmat(ibmatsub2,(iaconn(iatsub1)-1)*3+1)=bmn3x
         bmat(ibmatsub2,(iaconn(iatsub1)-1)*3+2)=bmn3y
         bmat(ibmatsub2,(iaconn(iatsub1)-1)*3+3)=bmn3z

      endif

      open(unit=20,file='bmat.dat',status='unknown')
      write(20,*)ncoord,3*natom
      do k=1,ncoord
         write(20,102)(bmat(k,j),j=1,3*natom)
      enddo
      close(20)

      if(ifile.eq.0)then
         open (unit=99,status='unknown')
         rewind (99)
         write (99,100) nameout
         rewind (99)
         read (99,101) command1
         close(99)
         call commrun(command1)         
      endif

 100  format ('cp -f bmat.dat output/'A8,'.dat')
 101  format (A100)
 102  format (10000(1X,1PE11.4))


cc now we can compute the C matrix using central difference
cc Cijk=d2qi/dxj/dxk

      jind=0
      kind=0
      jaind=0
      kaind=0

      do j=1,natomt
         if(idummy(j).ne.1)then
 
c*******************************************************************
cc perturb xj
            jind=jind+1
            jaind=jaind+1
            kind=0
            kaind=0
            do k=1,natomt
            if(idummy(k).ne.1)then
               kaind=kaind+1

cc  perturb xj + deltaxj
cc  perturb xk + deltaxk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               cooxpp(jaind)=cooxpp(jaind)+deltax
               cooxpp(kaind)=cooxpp(kaind)+deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj - deltaxj
cc  perturb xk + deltaxk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooxnp(jaind)=cooxnp(jaind)-deltax
               cooxnp(kaind)=cooxnp(kaind)+deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj + deltaxj
cc  perturb xk - deltaxk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               cooxpn(jaind)=cooxpn(jaind)+deltax
               cooxpn(kaind)=cooxpn(kaind)-deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj - deltaxj
cc  perturb xk - deltaxk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooxnn(jaind)=cooxnn(jaind)-deltax
               cooxnn(kaind)=cooxnn(kaind)-deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)

cc we can now compute the first of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltax/deltax
               enddo
c               stop
cc we now replicate for y displacement of k
cc  perturb xj + deltaxj
cc  perturb yk + deltayk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               cooxpp(jaind)=cooxpp(jaind)+deltax
               cooypp(kaind)=cooypp(kaind)+deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj - deltaxj
cc  perturb yk + deltayk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooxnp(jaind)=cooxnp(jaind)-deltax
               cooynp(kaind)=cooynp(kaind)+deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj + deltaxj
cc  perturb yk - deltayk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               cooxpn(jaind)=cooxpn(jaind)+deltax
               cooypn(kaind)=cooypn(kaind)-deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj - deltaxj
cc  perturb yk - deltayk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooxnn(jaind)=cooxnn(jaind)-deltax
               cooynn(kaind)=cooynn(kaind)-deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)

cc we can now compute the second of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltax/deltay
               enddo

cc we now replicate for z displacement of k
cc  perturb xj + deltaxj
cc  perturb zk + deltazk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               cooxpp(jaind)=cooxpp(jaind)+deltax
               coozpp(kaind)=coozpp(kaind)+deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj - deltaxj
cc  perturb zk + deltazk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooxnp(jaind)=cooxnp(jaind)-deltax
               cooznp(kaind)=cooznp(kaind)+deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj + deltaxj
cc  perturb zk - deltazk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               cooxpn(jaind)=cooxpn(jaind)+deltax
               coozpn(kaind)=coozpn(kaind)-deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb xj - deltaxj
cc  perturb zk - deltazk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooxnn(jaind)=cooxnn(jaind)-deltax
               cooznn(kaind)=cooznn(kaind)-deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc we can now compute the third of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltax/deltaz
               enddo
            endif
            enddo
c*******************************************************************
cc perturb yj
            jind=jind+1
            kind=0
            kaind=0
            do k=1,natomt
            if(idummy(k).ne.1)then
               kaind=kaind+1
cc  perturb yj + deltayj
cc  perturb xk + deltaxk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               cooypp(jaind)=cooypp(jaind)+deltay
               cooxpp(kaind)=cooxpp(kaind)+deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj - deltayj
cc  perturb xk + deltaxk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooynp(jaind)=cooynp(jaind)-deltay
               cooxnp(kaind)=cooxnp(kaind)+deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj + deltayj
cc  perturb xk - deltaxk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               cooypn(jaind)=cooypn(jaind)+deltay
               cooxpn(kaind)=cooxpn(kaind)-deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj - deltayj
cc  perturb xk - deltaxk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooynn(jaind)=cooynn(jaind)-deltay
               cooxnn(kaind)=cooxnn(kaind)-deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)

cc we can now compute the fourth of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltay/deltax
               enddo

cc we now replicate for y displacement of k
cc  perturb yj + deltayj
cc  perturb yk + deltayk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               cooypp(jaind)=cooypp(jaind)+deltay
               cooypp(kaind)=cooypp(kaind)+deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj - deltayj
cc  perturb yk + deltayk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooynp(jaind)=cooynp(jaind)-deltay
               cooynp(kaind)=cooynp(kaind)+deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj + deltayj
cc  perturb yk - deltayk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               cooypn(jaind)=cooypn(jaind)+deltay
               cooypn(kaind)=cooypn(kaind)-deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj - deltayj
cc  perturb yk - deltayk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooynn(jaind)=cooynn(jaind)-deltay
               cooynn(kaind)=cooynn(kaind)-deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc we can now compute the fifth of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltay/deltay
               enddo

cc we now replicate for z displacement of k
cc  perturb yj + deltayj
cc  perturb zk + deltazk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               cooypp(jaind)=cooypp(jaind)+deltay
               coozpp(kaind)=coozpp(kaind)+deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj - deltayj
cc  perturb zk + deltazk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooynp(jaind)=cooynp(jaind)-deltay
               cooznp(kaind)=cooznp(kaind)+deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj + deltayj
cc  perturb zk - deltazk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               cooypn(jaind)=cooypn(jaind)+deltay
               coozpn(kaind)=coozpn(kaind)-deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb yj - deltayj
cc  perturb zk - deltazk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooynn(jaind)=cooynn(jaind)-deltay
               cooznn(kaind)=cooznn(kaind)-deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc we can now compute the sixth of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltay/deltaz
               enddo
            endif
            enddo

c*******************************************************************
cc perturb zj
            jind=jind+1
            kind=0
            kaind=0
            do k=1,natomt
            if(idummy(k).ne.1)then
               kaind=kaind+1
cc  perturb zj + deltazj
cc  perturb xk + deltaxk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               coozpp(jaind)=coozpp(jaind)+deltaz
               cooxpp(kaind)=cooxpp(kaind)+deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj - deltazj
cc  perturb xk + deltaxk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooznp(jaind)=cooznp(jaind)-deltaz
               cooxnp(kaind)=cooxnp(kaind)+deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj + deltazj
cc  perturb xk - deltaxk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               coozpn(jaind)=coozpn(jaind)+deltaz
               cooxpn(kaind)=cooxpn(kaind)-deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj - deltazj
cc  perturb xk - deltaxk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooznn(jaind)=cooznn(jaind)-deltaz
               cooxnn(kaind)=cooxnn(kaind)-deltax
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)

cc we can now compute the seventh of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
c                  if(jind.eq.kind)then
c                     cmat(i,jind,kind)=(xintpp(i)-2*xint(i)+
c     $                    xintnn(i))/4./deltaz/deltax
c                  else
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltaz/deltax
c                  endif
               enddo

cc we now replicate for y displacement of k
cc  perturb zj + deltazj
cc  perturb yk + deltayk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               coozpp(jaind)=coozpp(jaind)+deltaz
               cooypp(kaind)=cooypp(kaind)+deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj - deltazj
cc  perturb yk + deltayk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooznp(jaind)=cooznp(jaind)-deltaz
               cooynp(kaind)=cooynp(kaind)+deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj + deltazj
cc  perturb yk - deltayk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               coozpn(jaind)=coozpn(jaind)+deltaz
               cooypn(kaind)=cooypn(kaind)-deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj - deltazj
cc  perturb yk - deltayk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooznn(jaind)=cooznn(jaind)-deltaz
               cooynn(kaind)=cooynn(kaind)-deltay
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc we can now compute the eigth of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltaz/deltay
c                  endif
               enddo

cc we now replicate for z displacement of k
cc  perturb zj + deltazj
cc  perturb zk + deltazk
               do ij=1,ncoord
                  xintpp(ij)=0.
               enddo
               do ik=1,natom
                  cooxpp(ik)=coox(ik)
                  cooypp(ik)=cooy(ik)
                  coozpp(ik)=cooz(ik)
               enddo
               coozpp(jaind)=coozpp(jaind)+deltaz
               coozpp(kaind)=coozpp(kaind)+deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
     $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj - deltazj
cc  perturb zk + deltazk
               do ij=1,ncoord
                  xintnp(ij)=0.
               enddo
               do ik=1,natom
                  cooxnp(ik)=coox(ik)
                  cooynp(ik)=cooy(ik)
                  cooznp(ik)=cooz(ik)
               enddo
               cooznp(jaind)=cooznp(jaind)-deltaz
               cooznp(kaind)=cooznp(kaind)+deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnp,cooynp,
     $     cooznp,xintnp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj + deltazj
cc  perturb zk - deltazk
               do ij=1,ncoord
                  xintpn(ij)=0.
               enddo
               do ik=1,natom
                  cooxpn(ik)=coox(ik)
                  cooypn(ik)=cooy(ik)
                  coozpn(ik)=cooz(ik)
               enddo
               coozpn(jaind)=coozpn(jaind)+deltaz
               coozpn(kaind)=coozpn(kaind)-deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxpn,cooypn,
     $     coozpn,xintpn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc  perturb zj - deltazj
cc  perturb zk - deltazk
               do ij=1,ncoord
                  xintnn(ij)=0.
               enddo
               do ik=1,natom
                  cooxnn(ik)=coox(ik)
                  cooynn(ik)=cooy(ik)
                  cooznn(ik)=cooz(ik)
               enddo
               cooznn(jaind)=cooznn(jaind)-deltaz
               cooznn(kaind)=cooznn(kaind)-deltaz
               call update_zmat(natom,natomt,intcoor,bislab,ibconn,
     $     iaconn,idconn,bname,anname,dname,atname,cooxnn,cooynn,
     $     cooznn,xintnn,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
     $     dconnt,atomlabel,ifilu)
cc we can now compute the ninth of nine components of the Cmatrix
               kind=kind+1
               do i=1,ncoord
                  if(abs(xintpp(i)-xintnp(i)).gt.300)then
                     if(xintnp(i).lt.0)then
                        xintnp(i)=xintnp(i)+360.
                     else if(xintnp(i).gt.0)then
                        xintnp(i)=xintnp(i)-360.
                     endif
                     if(abs(xintpp(i)-xintnp(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintpp(i)-xintpn(i)).gt.300)then
                     if(xintpn(i).lt.0)then
                        xintpn(i)=xintpn(i)+360.
                     else if(xintpn(i).gt.0)then
                        xintpn(i)=xintpn(i)-360.
                     endif
                     if(abs(xintpp(i)-xintpn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                  if(abs(xintnp(i)-xintnn(i)).gt.300)then
                     if(xintnn(i).lt.0)then
                        xintnn(i)=xintnn(i)+360.
                     else if(xintnn(i).gt.0)then
                        xintnn(i)=xintnn(i)-360.
                     endif
                     if(abs(xintnp(i)-xintnn(i)).gt.300)then
                        write(*,*)'something did not work here'
                        write(*,*)'kind= ',kind
                        write(*,*)'jind= ',jind
                        write(*,*)'intocoor= ',i
                        stop
                     endif
                  endif
                     cmat(i,jind,kind)=(xintpp(i)-xintnp(i)-xintpn(i)+
     $                    xintnn(i))/4./deltaz/deltaz
               enddo
            endif
            enddo
         endif
         write(*,*)'jind is ',jind
      enddo
      write(*,*)'kind is ',kind
      write(*,*)'jind is ',jind

      do k=1,ncoord
         do j=1,3*natom
            do i=1,3*natom
               if(iangind(k).eq.1)cmat(k,j,i)=cmat(k,j,i)/180.*pigr
            enddo
         enddo
      enddo


      open(unit=20,file='cmat.dat',status='unknown')
      write(20,*)ncoord,3*natom
      do i=1,ncoord
      write(20,*)i
         do j=1,3*natom
            write(20,105)(cmat(i,j,k),k=1,3*natom)
         enddo
      write(20,*)
      enddo
      close(20)

      if(ifile.eq.0)then
         open (unit=99,status='unknown')
         rewind (99)
         write (99,103) nameoutc
         rewind (99)
         read (99,104) command1
         close(99)
         call commrun(command1)         
      endif
 103  format ('cp -f cmat.dat output/'A8,'.dat')
 104  format (A100)
 105  format (10000(1X,1PE11.4))

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine intfreq_ts

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      include 'data_estoktp.fi'
      include 'param_estoktp.fi'

c      dimension natomnumb(natommx)
      dimension igrouptot(nhindmx),ipivotA(nhindmx),ipivotB(nhindmx)
      dimension ngroup(nhindmx,natommx)


      LOGICAL leof,lsec,ltit
      CHARACTER*1000 line,string
      CHARACTER*160 sename,word,word2,word3,title,title1,
     $ word4,word5,word6,word7

      character*20 filename
      character*20 stoichname
      character*180 command1
      character*30 cjunk

      include 'filcomm.f'

      command1='cp -f output/tsgta_bm.dat ./bmat.dat'
      call commrun(command1)

      command1='cp -f output/tsgta_cm.dat ./cmat.dat'
      call commrun(command1)

      open (unit=25,file='./data/reac1.dat',status='old')
      do while (WORD.NE.'NATOM')
         call LineRead (25)
         if (WORD.EQ.'END') then
            write (96,*) 'natom in reac1 must be defined'
            stop
         endif
      enddo
      read (25,*) natom1,natomt1
      close (25)

cc get data from react2 file
      if(iadd.eq.1.or.iabs.eq.1)then
         call LineRead (0)

         open (unit=25,file='./data/reac2.dat',status='old')

c     natomt is to account for dummy atoms
         do while (WORD.NE.'NATOM')
            call LineRead (25)
            if (WORD.EQ.'END') then
               write (96,*) 'natom must be defined'
               stop
            endif
         enddo
         read (25,*) natom2,natomt2
         close (25)
      endif

      if(iabs.eq.1.or.iadd.eq.1)then
         natom = natom1+natom2
      else if (iiso.eq.1.or.ibeta.eq.1) then
         natom = natom1
      endif

      if (iadd.eq.1) natomt = natomt1+natomt2
      if (iabs.eq.1) natomt = natomt1+natomt2+1
      if (iiso.eq.1) natomt = natomt1
      if (ibeta.eq.1) natomt = natomt1

      open (unit=15,file='./data/ts.dat',status='unknown')
      do while (WORD.NE.'NHIND')
         call LineRead (15)
         if (WORD.EQ.'END') then
            write (96,*) 'hind rotors must be defined'
            stop
         endif
      enddo
      read (15,*) nhind
      close (15)


      open (unit=133,file='RPHt_input_data1.dat',status='unknown')
c         open (unit=134,file='./data/hind_rot_head.dat',
c     +         status='unknown')

      write(133,*)'Number_of_Atoms: ',natom
      write(133,*)'Act_energy(kcal/mol):       0. '
      write(133,*)'Initial_Temperature:        200'
      write(133,*)'Temperature_steps:          40'
      write(133,*)'Temperature_increment:      40'
      write(133,*)'Delta_Energy_rea:           0.'
      write(133,*)'Delta_Energy_pro:           0.'
      write(133,*)'Maxstep:                    1'
      write(133,*)'Npointsint:                 5 '
      write(133,*)'Maxtdev:                    0.5'
      write(133,*)'Rearrange(1=yes,0=no)       1'
      write(133,*)'SaddlePoint                 1'
      write(133,*)'internalcoord(1=yes)     1'            
      write(133,*)'isct_vtst(1=vtst_sct,0=sct) 1'
      write(133,*)'zerocurvature(1)            0'
      write(133,*)'reduced_mass                1.0'
      write(133,*)'minimum_frequency            50'
      write(133,*)'anim_freq(if_Maxstep=1)       2'
      write(133,*)'onlyrotors(0=yes,1=no)        0'

      if(nhind.ne.0) then
         open (unit=15,file='./output/hrdata4proj_ts.dat'
     $        ,status='unknown')
         read (15,*)cjunk
         write (133,*)'proj_rea_coo(0=yes(def),1=no) ',1
         read (15,*)cjunk,nhind
         write (133,*)cjunk,nhind

         do ir=1,nhind
            read (15,*)cjunk,ipivotA(ir)
            write (133,*)cjunk,ipivotA(ir)
            read (15,*)cjunk,ipivotB(ir)
            write (133,*)cjunk,ipivotB(ir)
            read (15,*)cjunk,igrouptot(ir)
            write (133,*)cjunk,igrouptot(ir)
            read (15,*)cjunk,(ngroup(ir,igr),igr=1,
     $           igrouptot(ir))
            write (133,*)cjunk,(ngroup(ir,igr),igr=1,
     $           igrouptot(ir))

         enddo
         close(15)
      else
         if(inumpoints.eq.numpointsf+1)then
            write (133,*)'proj_rea_coo(0=yes(def),1=no) ',1
         else
            write (133,*)'proj_rea_coo(0=yes(def),1=no) ',iprojrcoo
         endif
         write (133,*)'numrotors ',0
      endif
      close(133)
