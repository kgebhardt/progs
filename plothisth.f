
      parameter(nmax=10000)
      real bl(nmax),bu(nmax),bin(nmax),bin2(nmax),binr(nmax)
      real xl(2),yl(2)

      sig=1.11

      xmin=-10.
      xmax=10.
      nbin=70
      xsize=(xmax-xmin)/float(nbin)
      do i=1,nbin
         bl(i)=xmin+float(i-1)*xsize
         bu(i)=bl(i)+xsize
         bin(i)=0.
      enddo

      open(unit=1,file='in',status='old')
      nt=0
      do i=1,100000000
         read(1,*,end=666) x1
         nt=nt+1
         do j=1,nbin
            if(x1.gt.bl(j).and.x1.le.bu(j)) bin(j)=bin(j)+1.
         enddo
      enddo
 666  continue
      close(1)

      ymin=1e10
      ymax=-ymin
      do i=1,nbin
         call getgausi(bl(i),bu(i),sig,gsum)
         gsum=gsum*float(nt)
         gsum=max(1.,gsum)
         bin2(i)=bin(i)-gsum
         binr(i)=bin2(i)/gsum
         if(bin(i).le.2.) binr(i)=0.
         print *,i,bin2(i),binr(i),bin(i),gsum
         ymin=min(ymin,binr(i))
         ymax=max(ymax,binr(i))
      enddo

      ymin=-0.1
      ymax=0.1

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.4)
      call pgslw(2)

      call pgenv(xmin,xmax,ymin,ymax,0,0)

      xl(1)=xmin
      xl(2)=xmax
      yl(1)=0.
      yl(2)=0.
      call pgline(2,xl,yl)

      call pgsci(1)
      do i=1,nbin
         if(binr(i).lt.0) call pgsci(2)
         xl(1)=bl(i)
         xl(2)=bl(i)
         yl(1)=0.
         yl(2)=binr(i)
         call pgline(2,xl,yl)
         xl(1)=bl(i)
         xl(2)=bu(i)
         yl(1)=binr(i)
         yl(2)=binr(i)
         call pgline(2,xl,yl)
         xl(1)=bu(i)
         xl(2)=bu(i)
         yl(1)=0.
         yl(2)=binr(i)
         call pgline(2,xl,yl)
         call pgsci(1)
      enddo

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
