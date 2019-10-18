
      parameter(nmax=10000)
      real x(nmax),y(nmax),ye(nmax),ya(nmax,100),xin(100)
      real yel(nmax),yeu(nmax),ydiff(nmax),ydiffn(nmax)
      character file1*80,file2*80,c1*3

c      call pgbegin(0,'?',3,3)
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xmin=1.
      xmax=1032.
      xmin=3500.
      xmax=5500.
      ymin=-1.5
      ymax=1.5
      ymin=0.
      ymax=0.25
      ymin=-40
      ymax=40.
c      xmin=1.
c      xmax=100.
      call pgsls(1)
      call pgslw(2)

      open(unit=1,file='in',status='old')

      nl=0
      ic=0
      do il=1,1000
c         do ia=1,12
c         do ia=1,70
         do ia=1,700
            read(1,*,end=666) file1
            open(unit=2,file=file1,status='old')
            if(ia.eq.1) then
               c1=file1(19:21)
               call pgsci(1)
               call pgenv(xmin,xmax,ymin,ymax,0,0)
c               call pglabel('Wavelength','Flux (1e-17 cgs)','')
c               call pglabel('Wavelength','Normalized Flux','')
               call pglabel('Wavelength','Offset from linear','')
c               call pglabel('Row','Counts','')
               call pgsch(2.0)
c               call pgmtxt('B',-1.4,0.5,0.5,c1)
               call pgsch(1.5)
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
            do i=1,n
               y(i)=y(i)-y(16)-20.
            enddo
            call pgslw(4)
            ic=ic+1
            if(ic.eq.13) ic=1
            call pgsci(ic)
            call pgline(n,x,y)
            call pgslw(1)
 555        continue
         enddo
      enddo
 666  continue
      close(1)

      call pgend

      end
