
      parameter(nmax=10000)
      real r(nmax),sb1(nmax),sb2(nmax),sb3(nmax)
      real rat2(nmax),rat3(nmax),rs(nmax),sig(nmax)
      real sig2(nmax),sig3(nmax)

      open(unit=1,file='out0',status='old')
      n=0
      xml=1.0
      fcon1=0.
c - tidal radii at 432 and 347
      fcon2=3.9
      fcon3=6.8
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3
         n=n+1
         r(n)=log10(x1)
         sb1(n)=log10((x3-fcon1)*xml)
         sb2(n)=log10(max(0.0001,(x3-fcon2)*xml))
         sb3(n)=log10(max(0.0001,(x3-fcon3)*xml))
c         rat2(n)=min(1.,1.-(x3-fcon2)/x3)
c         rat3(n)=min(1.,1.-(x3-fcon3)/x3)
         rat2(n)=max(0.,(x3-fcon2)/x3)
         rat3(n)=max(0.,(x3-fcon3)/x3)
c         print *,x1,sb1(n),sb2(n),sb3(n)
      enddo
 666  continue
      close(1)

      call pgbegin(0,'?',2,1)
      call pgpap(0.,0.5)
      call pgsch(1.4)
      call pgscf(2)
      call pgslw(2)

      xmin=log10(2.)
      xmax=log10(1000.)
      ymin=log10(0.9)
      ymax=log10(110.)
      call pgenv(xmin,xmax,ymin,ymax,0,30)
      call pglabel('R (arcsec)','\gS (L\D\(2281)\U/pc\U2\D)','')
      call pgsls(1)
      
      call pgsci(1)
      call pgline(n,r,sb1)
      call pgsci(2)
      call pgline(n,r,sb2)
      call pgsci(4)
      call pgline(n,r,sb3)
      call pgsci(1)

      open(unit=1,file='sigma.dat',status='old')
      ns=0
      sigc2=10.
      sigc3=11.4
      do i=1,nmax
         read(1,*,end=667) x1,x2,x3,x4
         ns=ns+1
         rs(ns)=x1
         sig(ns)=x2
         call xlinint(log10(x1),n,r,rat2,yp2)
         call xlinint(log10(x1),n,r,rat3,yp3)

         call xlinint(log10(x1),n,r,sb2,yp2s)
         call xlinint(log10(x1),n,r,sb3,yp3s)
         yp2s=10**yp2s
         signew=((yp2s+fcon2)*x2*x2-fcon2*sigc2*sigc2)/(yp2s+fcon2)
         signew=sqrt(max(0.,signew))
         sig2(ns)=signew

         yp3s=10**yp3s
         signew=((yp3s+fcon3)*x2*x2-fcon3*sigc3*sigc3)/(yp3s+fcon3)
         signew=sqrt(max(0.,signew))
         sig3(ns)=signew
         print *,x1,x2,sig2(ns),sig3(ns)
      enddo
 667  continue
      close(1)

      xmin=1.
      xmax=700.
      ymin=0.
      ymax=13.
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('R (arcsec)','Dispersion (km/s)','')
      call pgsls(1)
      
      call pgsci(1)
      call pgline(ns,rs,sig)
      call pgpt(ns,rs,sig,17)
      call pgsci(2)
      call pgline(ns,rs,sig2)
      call pgpt(ns,rs,sig2,17)
      call pgsci(4)
      call pgline(ns,rs,sig3)
      call pgpt(ns,rs,sig3,17)
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
      if(xp.gt.x(n)) yp=y(n)
      return
      end
