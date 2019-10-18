
      parameter(nmax=312*1032*112)
      parameter(nfmax=312)
      real arr1(1032,112),arr2(1032,112),xftf1(1032,112)
      real spec(1032*112*nfmax),wave(1032*112*nfmax)
      real wmaster(nmax),xmaster(nmax),waves(1036),xata(1036,312)
      real speca(1032,nfmax*112),wavea(1032,nfmax*112)
      real xftf(1032,nfmax*112),specs(1036,nfmax*112)
      real waveata(1036),xata1(1036),specsub(1032),wavesub(1032)
      real xin(1032*112),yin(1032*112),xp(10)
      real xin2(1032*112),yin2(1032*112),yin3(112)
      integer naxes(2),ifibc(112),ipeaks(10)
      character camp*2,cspecid*3,cifu*3,cifupos*3
      character camp2*2,cspecid2*3,cifu2*3,cifupos2*3
      character file1*100,cds9(312)*14,cmth*6

      iplot=0

      if(iplot.eq.1) then
         call pgbegin(0,'?',2,2)
c         call pgask(.false.)
         call pgpap(0.,1.)
         call pgscf(2)
         call pgsch(1.5)
         call pgslw(2)
      endif

c- set iuse to fit calibration from frame (0), or read in calibration (1)
      open(unit=1,file='vred2.in',status='old')
      read(1,*) imth,iuse
      close(1)
      write(cmth,2001) imth
 2001 format(i6)

      nwt=17
      wtmin=-1.0
      wtmax=1.0
      whalf=6.

c- read in all fits files
      open(unit=1,file='list',status='old')
      nt=0
      nta=0
      do iall=1,nfmax
         read(1,*,end=666) file1
         call getfits(file1,arr1,cspecid,camp,cifu,cifupos)
         write(file1(26:28),1001) "wav"
         call getfits(file1,arr2,cspecid2,camp2,cifu2,cifupos2)

c- read in amp-to-amp from twilights
         call readcal2(xftf1,cifupos,camp,waveata,xata1,iata,cmth,iuse)
         nta=nta+1
         cds9(nta)=cspecid//"_"//cifupos//"_"//cifu//"_"//camp
         do i=1,1032
            jin=1
            do j=1,112
               ntt=(nta-1)*112+j
               ftf0=xftf1(i,j)
               wavea(i,ntt)=arr2(i,j)
               if(iata.eq.1) arr1(i,j)=0.
               call xlinint2(wavea(i,ntt),1036,waveata,xata1,ata0,
     $              jin,jout)
               jin=max(1,jout-5)
               if(ftf0.eq.0) ftf0=1.
               if(ata0.eq.0) ata0=1.
               speca(i,ntt)=arr1(i,j)/ftf0/ata0
               if(speca(i,ntt).gt.-1e10.and.speca(i,ntt).lt.1e10) then
               else
                  speca(i,ntt)=0.
               endif
               if(arr1(i,j).ne.0) then
                  nt=nt+1
                  wave(nt)=wavea(i,ntt)
                  spec(nt)=speca(i,ntt)
               endif
            enddo
         enddo
      enddo
 666  continue
      close(1)
      print *,ntt,nta

c- make master sky
      call getmaster(nt,wave,spec,nmaster,wmaster,xmaster,1000,iuse,0)

c- find peaks for centering
      call getpeaks(nmaster,wmaster,xmaster,xp,ipeaks,iuse,0)
      call sort(4,xp)
      npeak=4
      xp(1)=4359.
      xp(2)=5084.
      xp(3)=5199.
      xp(4)=5461.

c- make new ftf and a2a
      call getftf(nmaster,wmaster,xmaster,nta,wavea,speca,xftf,
     $     waves,xata,iuse,0)
      if(iuse.eq.0) goto 555

c- iterate

c      nta=8
c- subtract master sky
      open(unit=12,file='amp.out',status='unknown')
      write(12,"(a45)") "Spc_slt_iid_am Factor N_c    Avg  Scale Woff"
      open(unit=13,file='ds9reg.dat',status='unknown')
      do iall=1,nta
c         print *,iall,nta
c- first get combined spectra for an amplifier
         nin=0
         do j=1,112
            ntt=(iall-1)*112+j
            do i=1,1032
               if(speca(i,ntt).ne.0.) then
                  nin=nin+1
                  xin(nin)=wavea(i,ntt)
                  yin(nin)=speca(i,ntt)
               endif
            enddo
         enddo
         call getmaster(nin,xin,yin,nin2,xin2,yin2,21,iuse,0)

c- find the continuum sources
         nin3=0
         do j=1,112
            ntt=(iall-1)*112+j
            nin=0
            jin=1
            do i=1,1032
               if(speca(i,ntt).ne.0.) then
                  nin=nin+1
                  nin3=nin3+1
                  call xlinint2(wavea(i,ntt),nin2,xin2,yin2,wv,jin,jout)
                  jin=jout
                  xin(nin3)=speca(i,ntt)-wv
                  yin(nin)=speca(i,ntt)-wv
               endif
            enddo
            call biwgt(yin,nin,xb,xs)
            yin3(j)=xb
         enddo
         call biwgt(xin,nin3,xb,xs)
         xcut=1.5*xs
         ncut=0
         do j=1,112
            ifibc(j)=0
            if(yin3(j).gt.xcut) ifibc(j)=1
            if(ifibc(j).eq.1) ncut=ncut+1
         enddo

         nin=0
         do j=1,112
            if(ifibc(j).eq.0) then
               ntt=(iall-1)*112+j
               do i=1,1032
                  if(speca(i,ntt).ne.0.) then
                     nin=nin+1
                     xin(nin)=wavea(i,ntt)
                     yin(nin)=speca(i,ntt)
                  endif
               enddo
            endif
         enddo
         call getmaster(nin,xin,yin,nin2,xin2,yin2,21,iuse,0)

         jin=1
         do i=1,nin2
            call xlinint2(xin2(i),nmaster,wmaster,xmaster,xv,jin,jout)
            jin=jout
            xin(i)=yin2(i)/xv
         enddo
         call biwgt(xin,nin2,xb,xs)
         if(ncut.ge.10) xb=1.
         if(xb.le.0.1) xb=1.

c- find the wavelength offset for the 4 lines
         xmin=1e10
         do iw=1,nwt
            wtry0=wtmin+float(iw-1)*(wtmax-wtmin)/float(nwt-1)
            nin=0
            do j=1,112
               ntt=(iall-1)*112+j
               jin=1
               if(ifibc(j).eq.0) then
                  do i=1,1032
                     wtry=wavea(i,ntt)+wtry0
                     stry=speca(i,ntt)
                     if(stry.ne.0.) then
                        do ip=1,npeak
                           if(wtry.gt.(xp(ip)-whalf).and.
     $                          wtry.lt.(xp(ip)+whalf)) then
                              call xlinint2(wtry,nmaster,wmaster,
     $                             xmaster,xv,jin,jout)
                              jin=jout
                              nin=nin+1
                              xin(nin)=(stry-xv*xb)**2
                           endif
                        enddo
                     endif
                  enddo
               endif
            enddo
            call biwgt(xin,nin,xbw,xsw)
            if(xbw.lt.xmin) then
               xmin=xbw
               wtmina=wtry0
            endif
         enddo

c- now subtract
         nin=0
         do j=1,112
            ntt=(iall-1)*112+j
            jin=1
            do i=1,1032
               wtry=wavea(i,ntt)+wtmina
               call xlinint2(wtry,nmaster,wmaster,xmaster,
     $              xv,jin,jout)
               wavesub(i)=wtry
               if(speca(i,ntt).eq.0.) then
                  specsub(i)=0.
               else
                  specsub(i)=speca(i,ntt)-xv*xb
               endif
               jin=jout
               if(ifibc(j).eq.0.and.specsub(i).ne.0.) then
                  nin=nin+1
                  xin(nin)=specsub(i)
               endif
            enddo
            jin=1
            do i=1,1036
               call xlinint2(waves(i),1032,wavesub,specsub,xv,jin,jout)
               specs(i,ntt)=xv
               jin=jout
            enddo
         enddo
         call biwgt(xin,nin,xbr,xsr)

c- write out the amplifier info
         write(12,1201) cds9(iall),xb,ncut,xbr,xsr,wtmina
         ix=10
         iy=iall*112-56
         ibad=0
         if(xb.lt.0.7) ibad=1
         if(xb.gt.1.1) ibad=1
         if(xsr.gt.20.) ibad=1
         if(xsr.lt.1.) ibad=1
         if(xsr.gt.16.and.abs(wtmina).gt.0.9) ibad=1
         if(ibad.eq.0) write(13,1301) ix,iy,cds9(iall)
         if(ibad.eq.1) write(13,1302) ix,iy,cds9(iall)
      enddo
      close(12)
      close(13)

 555  continue
c- write out each a2a from twilight
      if(iuse.eq.0) then
         open(unit=1,file='list',status='old')
         do iall=1,nta
            read(1,*,end=666) file1
            write(file1(26:32),1002) "ata.dat"
            open(unit=2,file=file1(1:32),status='unknown')
            do i=1,1036
               write(2,*) waves(i),xata(i,iall)
            enddo
            close(2)
         enddo
         close(1)
      endif

c- write out the sky-subtracted and rectified fits file
      if(iuse.eq.1) then
         naxis=2
         naxes(1)=1036
         naxes(2)=nta*112
         ier=0
         iblock=1
         igc=0
         call ftinit(51,'out.fits',iblock,ier)
         call ftphps(51,-32,naxis,naxes,ier)
         call ftp2de(51,igc,1036,naxes(1),naxes(2),specs,ier)
         call ftclos(51,ier)
      endif

 1001 format(a3)
 1002 format(a7)
 1201 format(a14,1x,f6.3,1x,i3,3(1x,f7.2))
 1301 format("# text("i2,", ",i5,") textangle=90 text={"a14"}")
 1302 format("# text("i2,", ",i5,
     $     ") textangle=90 text={"a14"} color=red")

      end

      subroutine getftf(nmaster,wmaster,xmaster,nta,wavea,speca,
     $     xftf,waves,xata,iuse,iplot)
      parameter(nmax=312*1032*112)
      parameter(nfmax=312)
      real wmaster(nmax),xmaster(nmax)
      real speca(1032,nfmax*112),wavea(1032,nfmax*112)
      real xftf(1032,nfmax*112),xin(nmax),yin(nmax)
      real yin2(nmax),yina(1032),ftfa(1036,112)
      real waves(1036),wave1(1032),spec1(1032),xata(1036,312)
      integer iftf(112)

      do i=1,1036
         waves(i)=3470.+float(i-1)*2.
      enddo

      if(iuse.eq.1) return

      do iall=1,nta
         if(iplot.eq.1) call pgenv(3470.,5540.,0.7,1.4,0,0)
         ic=0
         do j=1,112
            jt=(iall-1)*112+j
            iftf(j)=0
            nin=0
            do i=1,1032
               wave1(i)=wavea(i,jt)
               spec1(i)=speca(i,jt)
            enddo

c- first find which fibers to use
            jin=1
            jin2=1
            do i=1,1036
               call xlinint2b(waves(i),nmaster,wmaster,xmaster,xm1,
     $              jin,jout)
               jin=jout
               call xlinint2b(waves(i),1032,wave1,spec1,xs1,jin2,jout2)
               jin2=jout2
               ftfa(i,j)=xs1/xm1
               if(xs1.ne.0) then
                  nin=nin+1
                  xin(nin)=waves(i)
                  yin(nin)=xs1/xm1
               endif
            enddo
            do i=1,nin
               yin2(i)=yin(i)
            enddo
            call biwgt(yin2,nin,xb,xs)
c            print *,j,xb,xs
            if(xb.gt.0.65.and.xb.lt.1.3) then
               iftf(j)=1
               ic=ic+1
               if(ic.eq.14) ic=0
            else
               print *,j,xb,xs
            endif
         enddo

c- then compile to measure a2a
         do i=1,1036
            nin=0
            do j=1,112
               if(iftf(j).eq.1) then
                  nin=nin+1
                  yin(nin)=ftfa(i,j)
               endif
            enddo
            call biwgt(yin,nin,xb,xs)
            yin2(i)=xb
         enddo

         if(iplot.eq.1) then
            call pgslw(5)
            call pgsci(1)
            call pgline(1036,waves,yin2)
         endif

         do i=1,1036
            xata(i,iall)=yin2(i)
         enddo
      enddo

      return
      end

      subroutine getmaster(n,wave,spec,nmaster,wmaster,xmaster,ism,
     $     iuse,iplot)
      parameter(nmax=312*1032*112)
      real wave(n),spec(n)
      real wmaster(nmax),xmaster(nmax)

      call sort2(n,wave,spec)
      call smooth3(n,wave,spec,nmaster,wmaster,xmaster,ism)
c      print *,n,nmaster,float(nmaster)/1032.

      xmin=3500.
      xmax=5500.
      xmin=4035.
      xmax=4060.

      if(iplot.eq.1) then
         ymin=1e10
         ymax=-ymin
         do i=1,nmaster
            if(wmaster(i).gt.xmin.and.wmaster(i).lt.xmax) then
               ymin=min(ymin,xmaster(i))
               ymax=max(ymax,xmaster(i))
            endif
         enddo
         call pgenv(xmin,xmax,ymin,ymax,0,0)
         call pgline(nmaster,wmaster,xmaster)
         call pgsci(1)
      endif

      return
      end

      subroutine getfits(file1,arr,cspecid,camp,cifu,cifupos)
      real arr(1032,112)
      integer naxes(2)
      character file1*100,cname*15
      character camp*2,cspecid*3,cifu*3,cifupos*3
      logical simple,extend,anyf

      im1=0
      ier=0
      iext=1
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im1,iext,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,1032,ncol,nrow,arr,anyf,ier)

      ier2=0
      call ftgkys(im1,"SPECID",cspecid,cname,ier2)
      call ftgkys(im1,"AMP",camp,cname,ier2)
      call ftgkys(im1,"IFUID",cifu,cname,ier2)
      call ftgkys(im1,"IFUSLOT",cifupos,cname,ier2)

      call ftclos(im1,ier)
 706  continue
      if(ier.ne.0) print *,"No file for: ",file1
      return
      end

      subroutine smooth3(nin,xin,yin,nin2,xin2,yin2,ns)
      parameter(nmax=312*1032*112)
      real xin(nin),yin(nin),xin2(nmax),yin2(nmax)
      real yin3(nmax)

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
            yin3(n)=yin(j)
         enddo
         call biwgt(yin3,n,xb,xs)
         sum1=sum1/float(n)
         sum2=sum2/float(n)
         nin2=nin2+1
         xin2(nin2)=sum1
         yin2(nin2)=sum2
         yin2(nin2)=xb
      enddo

      return
      end

      subroutine xlinint(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.le.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.gt.x(n)) yp=y(n)
      return
      end

      subroutine xlinint2(xp,n,x,y,yp,jin,jout)
      real x(n),y(n)
c      do j=1,n-1
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

      subroutine xlinint2b(xp,n,x,y,yp,jin,jout)
      real x(n),y(n)
      do j=jin,n-1
         if(xp.ge.x(j).and.xp.le.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            jout=j
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.gt.x(n)) yp=y(n)
      jout=1
      return
      end

      subroutine readcal2(xftf,cifupos,camp,waveata,xata1,iata,
     $     cmth,iuse)

      real xftf(1032,112),waveata(1036),xata1(1036)
      integer naxes(2)
      character file1*100,file2*100
      character cifupos*3,camp*2,cmth*6
      logical simple,extend,anyf

c- get fiber-to-fiber
      file1="/work/03946/hetdex/maverick/virus_config/lib_calib/"//cmth
     $     //"/i"//cifupos//"a"//camp//"cmbf.fits"

      im1=0
      ier=0
      iext=1
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
         goto 707
      endif
      call ftmahd(im1,iext,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,1032,ncol,nrow,xftf,anyf,ier)
      call ftclos(im1,ier)

 707  continue
      if(ier.ne.0) print *,"No f2f for: ",file1

c- get amp-to-amp; if creating a2a (iuse=0) then return
      iata=0
      if(iuse.eq.0) return

      file1="/work/03946/hetdex/maverick/virus_config/lib_calib/"//cmth
     $     //"/i"//cifupos//"a"//camp//"ata.dat"

      open(unit=2,file=file1,status='old',err=708)
      do i=1,1036
         read(2,*) x1,x2
         waveata(i)=x1
         xata1(i)=x2
      enddo
      close(2)
      goto 709
 708  continue
      print *,"No a2a for: ",file1
      iata=1
      do i=1,1036
         waveata(i)=3470.+float(i-1)*2.
         xata1(i)=1.
      enddo
 709  continue

      return
      end

      subroutine getpeaks(nm,wm,xm,xp,ipeaks,iuse,iplot)
      parameter (narrm=50000)
      real wm(nm),xm(nm),xp(10)
      real xin(narrm),yin(narrm),yin2(narrm)
      integer ipeaks(10)

      if(iuse.eq.0) return

c- find peaks in master                                                                                                                                                     
      npeak=4
      call smooth3(nm,wm,xm,nin,xin,yin,401)
      ymin=1e10
      ymax=-ymin
      do i=1,nm
         call xlinint(wm(i),nin,xin,yin,yv)
         yin2(i)=xm(i)-yv
         ymin=min(ymin,yin2(i))
         ymax=max(ymax,yin2(i))
      enddo

      if(iplot.eq.1) then
         call pgenv(wm(1),wm(nm),ymin,ymax,0,0)
c         call pgenv(wm(1),wm(nm),50.,300.,0,0)
         call pgline(nm,wm,yin2)
c         call pgline(nm,wm,xm)
c         call pgsci(2)
c         call pgline(nin,xin,yin)
c         call pgsci(1)
      endif

      do ip=1,npeak
         ymax=-1e10
         do i=1,nm
            if(ip.gt.1) then
               do ipc=1,ip-1
                  if(abs(xp(ipc)-wm(i)).lt.15) goto 666
               enddo
            endif
            if(yin2(i).gt.ymax) then
               ymax=yin2(i)
               xmax=wm(i)
               imax=i
            endif
 666        continue
         enddo
         xp(ip)=xmax
         ipeaks(ip)=imax
      enddo
      
      return
      end
