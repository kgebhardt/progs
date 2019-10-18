
      parameter(nmax=10000)
      real bl(nmax),bu(nmax),bin(nmax),bin2(nmax),binr(nmax)
      real xl(2),yl(2),bing(nmax),bc(nmax)
      character file1*40

      sig=1.11

      xmin=-10.
      xmax=10.
      nbin=60
      xsize=(xmax-xmin)/float(nbin)
      do i=1,nbin
         bl(i)=xmin+float(i-1)*xsize
         bu(i)=bl(i)+xsize
         bc(i)=(bu(i)+bl(i))/2.
         bin(i)=0.
      enddo

      call pgbegin(0,'?',3,3)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.4)
      call pgslw(2)

      open(unit=2,file='listin',status='old')
      do ilist=1,1000
         read(2,*,end=667) file1

      open(unit=1,file='in',status='old')
      open(unit=1,file=file1,status='old')

      do i=1,nbin
         bin(i)=0.
      enddo
      nt=0
      do i=1,100000000
         read(1,*,end=666) x1,i2
         if(i2.lt.10) then
            nt=nt+1
            do j=1,nbin
               if(x1.gt.bl(j).and.x1.le.bu(j)) bin(j)=bin(j)+1.
            enddo
         endif
      enddo
 666  continue
      close(1)

      ymin=1e10
      ymax=-ymin
      do i=1,nbin
         call getgausi(bl(i),bu(i),sig,gsum)
         gsum=gsum*float(nt)
         bing(i)=gsum
c         print *,i,bin(i),bing(i)
         ymax=max(ymax,bin(i),bing(i))
         bin(i)=log10(max(0.1,bin(i)))
         bing(i)=log10(max(0.1,bing(i)))
      enddo

      ymin=log10(1.)
      ymax=log10(ymax)+0.1

      call pgsch(1.8)
      call pgenv(xmin,xmax,ymin,ymax,0,20)
      call pglabel(file1(1:15),"","")

c      call pgsci(15)
      call pgsci(1)
      call pgslw(1)
      do i=1,nbin
         xl(1)=bl(i)
         xl(2)=bl(i)
         yl(1)=0.
         yl(2)=bin(i)
         call pgline(2,xl,yl)
         xl(1)=bl(i)
         xl(2)=bu(i)
         yl(1)=bin(i)
         yl(2)=bin(i)
         call pgline(2,xl,yl)
         xl(1)=bu(i)
         xl(2)=bu(i)
         yl(1)=0.
         yl(2)=bin(i)
         call pgline(2,xl,yl)
      enddo

      call pgslw(5)
      call pgsci(2)
      call pgline(nbin,bc,bing)
      call pgslw(2)
      call pgsci(1)

      enddo
 667  continue
      close(2)

      call pgend

      end

      subroutine getgausi(b1,b2,sig,gsum)
      parameter(pi=3.14159)
      
      ng=100
      gsum=0.
      do i=1,ng
         r=b1+float(i-1)*(b2-b1)/float(ng-1)
         r=(r/sig)**2
         gsum=gsum+exp(-r/2.)
      enddo
      gsum=gsum*(b2-b1)/float(ng-1)/sqrt(2.*pi)/sig

      return
      end
