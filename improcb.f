
      parameter (narrm=2000)
      real xd(narrm,narrm),arr(narrm,narrm),xin(narrm*narrm)
      real xov(narrm),xbias(narrm,narrm)
      integer naxes(2)
      character file1*40,cspec*100
      logical simple,extend,anyf

      file1='in.fits'
      iext=1

      im1=0
      ier=0
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im1,iext,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      if(naxes(1).gt.narrm.or.naxes(2).gt.narrm) then
         write(*,"('Arrays too small - make narrm bigger')")
         goto 706
      endif
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,narrm,ncol,nrow,xd,anyf,ier)

      ix1=1034
      ix2=1060
      iy1=10
      iy2=1032

      ixt1=1
      ixt2=1032
      iyt1=1
      iyt2=1032

      ixn1=300
      ixn2=700
      iyn1=300
      iyn2=700

      n=0
      nov=0
      do j=iy1,iy2
         nov=nov+1
         do i=ix1,ix2
            n=n+1
            xin(n)=xd(i,j)
         enddo
      enddo
      call biwgt(xin,n,xb,xs)
      xov(nov)=xb
      xbo=xb

      ny=0
      do j=iyt1,iyt2
         ny=ny+1
         nx=0
         do i=ixt1,ixt2
            nx=nx+1
            arr(nx,ny)=xd(i,j)-xb
         enddo
      enddo      

      naxes(1)=ixt2-ixt1+1
      naxes(2)=iyt2-iyt1+1

      n=0
      do j=iyn1,iyn2
         do i=ixn1,ixn2
            n=n+1
            xin(n)=arr(i,j)
         enddo
      enddo      

      call biwgt(xin,n,xb,xs)
      print *,'Overscan, Average count is : ',xbo,xb,xs
      open(unit=11,file='out',status='unknown')
      write(11,*) xbo,xb,xs
      close(11)

      ier=0
      call ftclos(im1,ier)
      call ftinit(51,'improc.fits',iblock,ier)
      call ftphps(51,-32,naxis,naxes,ier)
      if(ier.ne.0) then
         print *,'Error in output file ',ier
         goto 706
      endif
      call ftp2de(51,igc,narrm,naxes(1),naxes(2),arr,ier)
      call ftclos(51,ier)

 706  continue
      end

