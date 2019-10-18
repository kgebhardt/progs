
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax)
      real xb(nmax),yb(nmax)
      character file1*80,file2*80,c1*18

c      ibin=13
      ibin=1
      ib1=(ibin-1)/2
      xib=float(ibin)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xmin=22000.
      xmax=24000.
      xmin=8400.
      xmax=8800.
c      xmin=1.
c      xmax=112.
c      ymin=-0.03
c      ymax=0.03
      ymin=0.
      ymax=0.18
      ymin=600.
      ymax=1150.
      ymin=-0.4
      ymax=0.4
      ymin=0.65
      ymax=1.05

      call pgsls(1)
      call pgslw(2)
      call pgsci(1)
      call pgsch(1.2)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
c      call pglabel('Wavelength','Residual from mastersky (600 cts)','')
c      call pglabel('Fiber','FtF correction','')
      call pgsch(1.4)
c      call pglabel('Wavelength','Throughput for 50 sq-m','')
c      call pglabel('Wavelength','Resolving Power','')
c      call pglabel('Fiber','Fiber Renormalization','')
c      call pglabel('Fiber','Wavelength Offset','')
c      call pglabel('Fiber','Sigma instrumental','')
      call pgslw(3)
      call pgsch(1.)

      open(unit=1,file='list',status='old')

      nl=0
      ic=0
      do il=1,10000
         read(1,*,end=666) file1
         open(unit=2,file=file1,status='old')
         n=0
         do i=1,nmax
            read(2,*,end=667) x1,x2
            if(x2.ge.-1e10.or.x2.le.1e10) then
               n=n+1
               x(n)=x1
               y(n)=x2
            else
               goto 866
            endif
c            if(x(n).eq.4740.) print *,x2,file1
         enddo
 667     continue
         close(2)
         nbb=0
         do j=1,n,ibin
            nbb=nbb+1
            istart=max(0,j-ib1)
            iend=istart+ibin-1
            if(iend.gt.n) then
               iend=n
               istart=n-ibin+1
            endif
            sum=0.
            nb=0
            do is=istart,iend
               sum=sum+y(is)
               nb=nb+1
               yb(nb)=y(is)
               xb(nb)=x(is)
            enddo
            call biwgt(yb,nb,xbb,xsb)
            yn(nbb)=xbb
            call biwgt(xb,nb,xbb,xsb)
            xn(nbb)=xbb
         enddo
         c1=file1(1:17)
         ic=ic+1
         icp=int(float(ic)/4.)+1
         if(icp.eq.14) ic=1
         if(ic.eq.14) ic=1
         call pgsci(icp)
         call pgsci(ic)
c         call pgline(n,x,y)
         call pgline(nbb,xn,yn)
 866     continue
         close(2)
      enddo
 666  continue
      close(1)

      call pgend

      end
