
      parameter (narrm1=2000,narrm2=2000)
      real xd(narrm1,narrm2),xin(narrm1),xrow(narrm1),xsing(narrm1)
      real xback(5),yback(5),xbacka(narrm1,narrm1),ybacka(narrm1,narrm1)
      integer naxes(2),ixl(narrm1),ixu(narrm1),ixm1(narrm1),ixm2(narrm1)
      character file1*120
      logical simple,extend,anyf

      file1='in.fits'
      iext1=1
      ismooth=13

      im1=0
      ier=0
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im1,iext1,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      if(naxis.eq.1) naxes(2)=1
      if(naxes(1).gt.narrm1.or.naxes(2).gt.narrm2) then
         write(*,"('Arrays too small - make narrm bigger')")
         goto 706
      endif
      ncol=naxes(1)
      nrow=max(1,naxes(2))
      call ftg2de(im1,igc,-666.,narrm1,ncol,nrow,xd,anyf,ier)
      call ftclos(im1,ier)

      nrcheck=35
c- first do bottom
      do i=1,ncol
         n=0
         do j=1,nrcheck
            n=n+1
            xin(n)=xd(i,j)
            if(n.gt.4) then
               call biwgt(xin,n,xb,xs)
               xstd=(xd(i,j)-xb)/xs
               if(xstd.gt.6.) then 
                  xrow(i)=float(j)
                  goto 555
               endif
            endif
         enddo
         if(i.gt.1) xrow(i)=xrow(i-1)
 555     continue
      enddo

      ibin=9
      ihalf=nint(float(ibin+1)/2.)
      do i=1,ncol
         ibs=max(1,i-ihalf)
         ibe=ibs+ibin-1
         if(ibe.gt.ncol) then
            ibs=ncol-ibin+1
            ibe=ncol
         endif
         n=0
         do j=ibs,ibe
            n=n+1
            xin(n)=xrow(j)
         enddo
         call biwgt(xin,n,xb,xs)
         ixl(i)=nint(xb-1)
      enddo

c- now the top
      do i=1,ncol
         n=0
         do j=nrow,nrow-nrcheck,-1
            n=n+1
            xin(n)=xd(i,j)
            if(n.gt.4) then
               call biwgt(xin,n,xb,xs)
               xstd=(xd(i,j)-xb)/xs
               if(xstd.gt.6.) then 
                  xrow(i)=float(j)
                  goto 556
               endif
            endif
         enddo
c         if(i.gt.1) xrow(i)=xrow(i-1)
         if(i.gt.1) xrow(i)=nrow
 556     continue
      enddo

      ibin=9
      ihalf=nint(float(ibin+1)/2.)
      do i=1,ncol
         ibs=max(1,i-ihalf)
         ibe=ibs+ibin-1
         if(ibe.gt.ncol) then
            ibs=ncol-ibin+1
            ibe=ncol
         endif
         n=0
         do j=ibs,ibe
            n=n+1
            xin(n)=xrow(j)
         enddo
         call biwgt(xin,n,xb,xs)
         ixu(i)=nint(xb+1)
      enddo

c- now the gaps in middle
      nhalf=ncol/2
      do i=1,ncol
         smin=1e10
         do j=ixl(i)+40,nhalf
            irs=j-4
            ire=j+4
            sum=0.
            do k=irs,ire
               sum=sum+xd(i,k)
            enddo
            if(sum.lt.smin) then
               smin=sum
               jmin=j
            endif
         enddo
         ixm1(i)=jmin
         smin2=1e10
         do j=nhalf,ixu(i)-40
            irs=j-4
            ire=j+4
            sum=0.
            do k=irs,ire
               sum=sum+xd(i,k)
            enddo
            if(sum.lt.smin2) then
               smin2=sum
               jmin=j
            endif
         enddo
         ixm2(i)=jmin
      enddo

c- now make the background array
      do i=1,ncol
         n=0
         js=1
         if(ixl(i).ge.11) js=6
         do j=js,ixl(i)
            n=n+1
            xin(n)=xd(i,j)
         enddo
         call biwgt(xin,n,xb,xs)
         xback(1)=xb
         yback(1)=float(ixl(i)+1)/2.
         if(n.lt.5) then
            xback(1)=0.
            yback(1)=1.
         endif
         n=0
         do j=ixm1(i)-4,ixm1(i)+4
            n=n+1
            xin(n)=xd(i,j)
         enddo
         call biwgt(xin,n,xb,xs)
         xback(2)=xb
         yback(2)=float(ixm1(i))
         n=0
         do j=ixm2(i)-4,ixm2(i)+4
            n=n+1
            xin(n)=xd(i,j)
         enddo
         call biwgt(xin,n,xb,xs)
         xback(3)=xb
         yback(3)=float(ixm2(i))
         n=0
         do j=ixu(i),nrow
            n=n+1
            xin(n)=xd(i,j)
         enddo
         call biwgt(xin,n,xb,xs)
         xback(4)=xb
         yback(4)=float(ixu(i)+nrow)/2.
         if(n.lt.5) then
            xback(4)=0.
            yback(4)=float(nrow)
         endif

c         if(xback(1).gt.6) xback(1)=0.
c         if(xback(4).gt.6) xback(4)=0.

         if(xback(1).eq.0) xback(1)=xback(2)
         if(xback(4).eq.0) xback(4)=xback(3)

         frac=0.1
         ybnew=yback(1)+frac*(yback(2)-yback(1))
         xbnew=(1.-frac)*xback(2)+frac*xback(1)
         xback(5)=xback(4)
         yback(5)=yback(4)
         xback(4)=xback(3)
         yback(4)=yback(3)
         xback(3)=xback(2)
         yback(3)=yback(2)
         xback(2)=xbnew
         yback(2)=ybnew         
c         print *,i,xback
c         print *,yback

         do j=1,5
            ybacka(i,j)=yback(j)
            xbacka(i,j)=xback(j)
         enddo
      enddo

c - do the bottom rows individually

      isub=3
      if(ixl(500).gt.20) then
         isingle=17
      else
         isingle=max(2,ixl(500)-isub)
      endif
      do j=1,isingle
         nin=0
         do i=400,600
            nin=nin+1
            xin(nin)=xd(i,j)
         enddo
         call biwgt(xin,nin,xb,xs)
         if(xb.lt.-10.) xb=0.
         xsing(j)=xb
c         print *,j,xb
      enddo   
      do i=1,ncol
         ybacka(i,2)=float(isingle)
         xbacka(i,2)=xsing(isingle)
      enddo

      do j=3,5
         do i=1,ncol
            xin(i)=xbacka(i,j)
         enddo
         call biwgt(xin,ncol,xb,xs)
         if(xb.gt.10) then
            do i=1,ncol
               xbacka(i,j)=xbacka(i,j-1)
            enddo
         endif
      enddo
      do j=1,5
c         print *,j,xbacka(500,j)
      enddo

      do j=1,5
         do i=1,ncol
            xin(i)=xbacka(i,j)
         enddo
         call biwgt(xin,ncol,xb,xs)
         do i=1,ncol
            xin(i)=xbacka(i,j)
            sd=(xin(i)-xb)/xs
            if(sd.gt.5.) then
               xin(i)=-666
            endif            
         enddo
         call smxin(ncol,xin,ismooth)
         do i=1,ncol
            xbacka(i,j)=xin(i)
c            print *,xin(i)
         enddo
      enddo

      do i=1,ncol
         do j=1,5
            xback(j)=xbacka(i,j)
            yback(j)=ybacka(i,j)
         enddo
         do j=1,nrow
            call xlinint2(float(j),5,yback,xback,xout)
            xd(i,j)=xout
         enddo
         do j=1,isingle
            xd(i,j)=xsing(j)
         enddo
      enddo

      naxis=2
      naxes(1)=ncol
      naxes(2)=nrow
      iblock=1
      igc=0
      ier=0
      call ftinit(50,'image.fits',iblock,ier)
      call ftphps(50,-32,naxis,naxes,ier)
      if(ier.ne.0) then
         print *,'Error in output file ',ier
      endif
      call ftp2de(50,igc,narrm1,naxes(1),naxes(2),xd,ier)
      call ftclos(50,ier)

 706  continue
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

      subroutine xlinint2(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            if(j.eq.2) then
               fac=(xp-x(j))/(x(j+1)-x(j))
               fac=fac**0.3
               yp=y(j)+(y(j+1)-y(j))*fac
            endif
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.gt.x(n)) yp=y(n)
      return
      end

      subroutine smxin(n,x,ibin)
      parameter(nmax=10000)
      real x(n),x2(nmax),xi2(nmax)
      real xb(nmax),yb(nmax),xn(nmax),yn(nmax)
      n2=0
      do i=1,n
         if(x(i).ne.-666) then
            n2=n2+1
            x2(n2)=x(i)
            xi2(n2)=float(i)
         endif
      enddo

      ib1=(ibin-1)/2
      xib=float(ibin)
      nbb=0
      do j=1,n2,ibin
         nbb=nbb+1
         istart=max(1,j-ib1)
         iend=istart+ibin-1
         if(iend.gt.n2) then
            iend=n2
            istart=n2-ibin+1
         endif
         sum=0.
         nb=0
         do is=istart,iend
            sum=sum+x2(is)
            nb=nb+1
            yb(nb)=x2(is)
            xb(nb)=xi2(is)
         enddo
         call biwgt(yb,nb,xbb,xsb)
         yn(nbb)=xbb
         call biwgt(xb,nb,xbb,xsb)
         xn(nbb)=xbb
c         print *,nbb,xn(nbb),yn(nbb)
      enddo

      do i=1,n
         call xlinint(float(i),nbb,xn,yn,xbv)
         x(i)=xbv
      enddo

      return
      end
