
      parameter(narrm=1680,nfmax=180)
      real xda(narrm,narrm,nfmax),xd(narrm,narrm)
      real xds(narrm,narrm),xin(narrm)
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

      open(unit=1,file='list',status='old')

      ntot=201
      n=0
      do ia=1,ntot
         read(1,*,end=666) file1
         n=n+1
         if(n.eq.ntot) then
            print *,"Too many files. Use imcmbbs_slow or increase array"
            goto 706
         endif

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
         call ftg2de(im1,igc,0.,narrm,ncol,nrow,xd,anyf,ier)
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
            do ia=1,n
               xin(ia)=xda(i,j,ia)
            enddo
            call biwgt(xin,n,xb,xs)
            xd(i,j)=xb
            xds(i,j)=xs
         enddo
      enddo

      call ftinit(51,'out1.fits',iblock,ier)
      call ftphps(51,-32,naxis,naxes,ier)
      call ftp2de(51,igc,narrm,naxes(1),naxes(2),xd,ier)
      call ftclos(51,ier)

      call ftinit(51,'out2.fits',iblock,ier)
      call ftphps(51,-32,naxis,naxes,ier)
      call ftp2de(51,igc,narrm,naxes(1),naxes(2),xds,ier)
      call ftclos(51,ier)

 706  continue
      end
