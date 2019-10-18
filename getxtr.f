
      real xda(1032,112,100),xd(1032,112)
      real xda2(1032,112,100),xds(1032,112)
      real xin(10000),yin(10000)
      real xin2(2000),yin2(2000),yin3(2000)
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

      call pgbegin(0,'?',2,2)
c      call pgask(.false.)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      open(unit=1,file='list',status='old')

      n=0
      do ia=1,100
         read(1,*,end=666) file1
         n=n+1

         im1=0
         ier=0
         iread=0
         iext=1
         call ftgiou(im1,ier)
         call ftopen(im1,file1,iread,iblock,ier)
         call ftmahd(im1,iext,ihd,ier)
         call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
         ncol=naxes(1)
         nrow=naxes(2)
         call ftg2de(im1,igc,0.,1032,ncol,nrow,xd,anyf,ier)
         call ftclos(im1,ier)
         do i=1,ncol
            do j=1,nrow
               xda(i,j,n)=xd(i,j)
            enddo
         enddo
      enddo
 666  continue

      do i=1,ncol
         do j=1,nrow
            xd(i,j)=xda(i,j,1)
         enddo
      enddo

      ncbin=15
      xbin=float(1032-1)/float(ncbin)
      ibin=nint(xbin)
      ibinh=nint(xbin/2.)
      istart=30
      do ia=1,n
         call pgslw(2)
c         call pgenv(1.,112.,-0.2,0.1,0,0)
         call pgenv(1.,112.,-5.,5.,0,0)

         do icol=istart,1032,ibin
            ilo=max(1,icol-ibinh)
            ihi=min(1032,icol+ibinh)
            if(ihi.eq.1030) ihi=1032
            icen=nint(float(ilo)+float(ihi)/2.)
            nin=0
            do j=1,nrow
               nin=nin+1
               xin(nin)=float(j)
               sum=0.
               nsum=0
               do isum=ilo,ihi
                  sum=sum+xda(isum,j,ia)-xd(isum,j)
                  nsum=nsum+1
               enddo
               yin(nin)=sum/float(nsum)
            enddo

c            call pgpt(nin,xin,yin,17)
            call smooth2(nin,xin,yin,nin2,xin2,yin2,9)
            call smooth(nin2,xin2,yin2,nin,xin,yin,yin3,0.)
            call pgslw(4)
            call pgsci(2)
c            call pgline(nin2,xin2,yin2)
            call pgsci(4)
            call pgline(nin,xin,yin)
            call pgsci(1)

            jin=1
            do j=1,nrow
c               call xlinint(float(j),nin,xin,yin,yv)
               call xlinint2(float(j),nin,xin,yin,yv,jin,jout)
               jin=jout
               do isum=ilo,ihi
                  xda2(isum,j,ia)=xda(isum,j,ia)-yv
               enddo
            enddo
         enddo
      enddo

      ncbin=15
      xbin=float(1032-1)/float(ncbin)
      ibin=nint(xbin)
      ibinh=nint(xbin/2.)
      istart=30
      do ia=1,n
         call pgslw(2)
         call pgenv(1.,112.,-0.2,0.1,0,0)
c         call pgenv(1.,112.,-5.,5.,0,0)

         do icol=istart,1032,ibin
            ilo=max(1,icol-ibinh)
            ihi=min(1032,icol+ibinh)
            icen=nint(float(ilo)+float(ihi)/2.)
            nin=0
            do j=1,nrow
               nin=nin+1
               xin(nin)=float(j)
               sum=0.
               nsum=0
               do isum=ilo,ihi
                  sum=sum+xda2(isum,j,ia)-xd(isum,j)
                  nsum=nsum+1
               enddo
               yin(nin)=sum/float(nsum)
            enddo

            call pgpt(nin,xin,yin,17)
         enddo
      enddo
      call pgend

      do i=1,1032
         do j=1,112
            do ia=1,n
               xin(ia)=xda2(i,j,ia)
            enddo
            call biwgt(xin,n,xb,xs)
            xd(i,j)=xb
            xds(i,j)=xs
         enddo
      enddo

      naxes(1)=1032
      naxes(2)=112
      call ftinit(51,'out1.fits',iblock,ier)
      call ftphps(51,-32,naxis,naxes,ier)
      call ftp2de(51,igc,1032,naxes(1),naxes(2),xd,ier)
      call ftclos(51,ier)

      call ftinit(51,'out2.fits',iblock,ier)
      call ftphps(51,-32,naxis,naxes,ier)
      call ftp2de(51,igc,1032,naxes(1),naxes(2),xds,ier)
      call ftclos(51,ier)

 706  continue
      end

      subroutine smooth2(nin,xin,yin,nin2,xin2,yin2,ns)
      parameter(narrm=2000)
      real xin(nin),yin(nin),xin2(narrm),yin2(narrm)

      nin2=0
      do i=1,nin,ns
         n=0
         sum1=0.
         sum2=0.
         jmin=i
         jmax=min(nin,i+ns)
         do j=jmin,jmax
            n=n+1
            sum1=sum1+xin(j)
            sum2=sum2+yin(j)
         enddo
         sum1=sum1/float(n)
         sum2=sum2/float(n)
         nin2=nin2+1
         xin2(nin2)=sum1
         yin2(nin2)=sum2
      enddo

      return
      end

      subroutine smooth(n,x,y,n2,x2,y2,y3,sm)
      parameter(nmax=20000,mm=2,nwk=nmax+6*(nmax*mm+1),mm2=mm*2)
      real x(n),y(n),x2(n2),y2(n2),y3(n)
      real*8 dx(nmax),dy(nmax),wx(nmax),cf(nmax),wk(nwk),val,splder
      real*8 q(mm2)

      if(n.gt.nmax) print *,'make nmax bigger in smooth'

      val=0.d0
      val=0.001d0
      val=dble(sm)
      md=3
      if(val.eq.0.) md=2
      m=2

      do i=1,n
         dx(i)=dble(x(i))
         dy(i)=dble(y(i))
         wx(i)=1.d0
      enddo

      call gcvspl(dx,dy,nmax,wx,1.d0,m,n,1,md,val,cf,nmax,wk,ier)
      if(ier.ne.0) print *,'ier= ',ier

      do i=1,n2
         in=i
         y2(i)=sngl(splder(0,m,n,dble(x2(i)),dx,cf,in,q))
      enddo

      do i=1,n
         in=i
         y3(i)=sngl(splder(0,m,n,dble(x(i)),dx,cf,in,q))
      enddo

      return
      end

      subroutine xlinint(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.ge.x(n)) yp=y(n)
      return
      end

      subroutine xlinint2(xp,n,x,y,yp,jin,jout)
      real x(n),y(n)
      do j=jin,n-1
         if(xp.ge.x(j).and.xp.le.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            jout=j
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=0.
      if(xp.gt.x(n)) yp=0.
      jout=1
      return
      end
