
      parameter(nmax=10000)
      real x(nmax),y(nmax),ye(nmax),ya(nmax,1036),xin(nmax)
      real yel(nmax),yeu(nmax),ydiff(nmax),ydiffn(nmax)
      real yin(nmax),yl(nmax),yh(nmax),ya2(nmax,1036)
      character file1*80,file2*80,c1*3

c      call pgbegin(0,'?',3,3)
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xmin=3500.
      xmax=5500.
      ymin=-15.
      ymax=15.
      call pgsls(1)
      call pgslw(2)

      open(unit=1,file='in',status='old')

      nl=0
      ic=0
      nall=0
      do il=1,1000
         do ia=1,700
            read(1,*,end=666) file1
            open(unit=2,file=file1,status='old')
            if(ia.eq.1) then
               c1=file1(19:21)
c               call pgsci(1)
c               call pgenv(xmin,xmax,ymin,ymax,0,0)
c               call pglabel('Wavelength','Offset from linear','')
c               call pgsch(1.5)
            endif
            n=0
            do i=1,8000
               read(2,*,end=667) x1,x2
               n=n+1
               x(n)=x1
               y(n)=x2
            enddo
 667        continue
            close(2)
            if(y(1).gt.10.) goto 555
            if(y(1).lt.-22.) goto 555
            if(y(416).gt.6.) goto 555
            if(y(1020).lt.-3.) goto 555
            y516=y(516)
            nall=nall+1
            do i=1,n
               y(i)=y(i)-y516-10.
               ya(nall,i)=y(i)
            enddo
            ic=ic+1
            if(ic.eq.13) ic=1
c            call pgslw(4)
c            call pgsci(ic)
c            call pgline(n,x,y)
c            call pgslw(1)
 555        continue
         enddo
      enddo
 666  continue
      close(1)

      do i=1,1036
         do ia=1,nall
            xin(ia)=ya(ia,i)
         enddo
         call biwgt(xin,nall,xb,xs)
         yin(i)=xb
         yl(i)=xb-xs
         yh(i)=xb+xs
c         print *,x(i),yin(i),xs
      enddo

      nw=1000
      ws=-2.
      we=2.
      do ia=1,nall
         rmsmin=1e10
         do iw=1,nw
            woff=ws+(we-ws)*float(iw-1)/float(nw-1)
            rms=0.
            do i=1,1036
               rms=rms+(yin(i)+woff-ya(ia,i))**2
            enddo
            if(rms.lt.rmsmin) then
               woffmin=woff
               rmsmin=rms
            endif
         enddo
         do i=1,1036
            ya2(ia,i)=ya(ia,i)-woffmin
         enddo
      enddo

      do i=1,1036
         do ia=1,nall
            xin(ia)=ya2(ia,i)
         enddo
         call biwgt(xin,nall,xb,xs)
         yin(i)=xb
         yl(i)=xb-xs
         yh(i)=xb+xs
         print *,x(i),yin(i),xs
      enddo

      call pgsci(1)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('Wavelength','Offset from linear','')
      call pgsci(1)
      call pgslw(5)
      call pgline(n,x,yin)
      call pgsls(4)
      call pgslw(3)
      call pgline(n,x,yl)
      call pgline(n,x,yh)

      call pgend

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
