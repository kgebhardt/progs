
      parameter (narrm=3000,nmax=1000)
      real xd(narrm,narrm),wave(narrm),xp(narrm),xerr2(narrm)
      real xdw(narrm,narrm),xdf(narrm,narrm),xdss(narrm,narrm)
      real xderr(narrm,narrm)
      real xsp(narrm),waved(narrm),wtp(narrm),tp(narrm)
      real xsky(narrm,narrm),sky(narrm),xspa(narrm,narrm)
      real xtrace(narrm,narrm),trace(narrm),xspab(narrm,narrm)
      real xflag(narrm,narrm),flag(narrm),x2dsub(narrm, narrm)
      real wa(narrm),fa(narrm),f2f(narrm),xerr(narrm),xin(100000)
      real wbadl(narrm),wbadu(narrm),xdout(1036,35000)
      real wbadla(narrm),wbadua(narrm)
      integer naxes(2),iwave(narrm)
      integer ixbs(1000),ixbe(1000),iybs(1000),iybe(1000)
      integer ixbsa(1000),ixbea(1000),iybsa(1000),iybea(1000)
      character file1*130,file2*130,file3*130,amp*14,a1*14,ae*5,aef*5
      character date*8,shot*3,abad(nmax)*14
      logical simple,extend,anyf

      nfib=112
      w0=4505.
      ww=1035.

      wmin=w0-ww
      wmax=w0+ww
      wsp=2.
      nw=nint((wmax-wmin)/wsp)+1
      do i=1,nw
         wave(i)=wmin+float(i-1)*wsp
      enddo

      convfac=(6.626e-27)*(3.e18)/360./5.e5

c - getting the throughput

c      open(unit=1,file=file2,status='old',err=576)
c      ntp=0
c      do i=1,narrm
c         read(1,*,end=676) x1,x2
c         ntp=ntp+1
c         wtp(ntp)=x1
c         tp(ntp)=x2
c      enddo
c 676  continue
c      close(1)
c      if(ntp.le.1) goto 576
c      goto 578
c 576  continue
c      close(1)

      open(unit=1,file=
     $     "/work/00115/gebhardt/maverick/detect/tp/tpavg.dat",
     $     status='old')
      ntp=0
      do i=1,narrm
         read(1,*,end=577) x1,x2
         ntp=ntp+1
         wtp(ntp)=x1
         tp(ntp)=x2
      enddo
 577  continue
      close(1)
 578  continue

c - get the bad pixel list
      open(unit=1,file='/work/03946/hetdex/hdr1/calib/badpix.new'
     $     ,status='old')
      nbada=0
      nbadwa=0
      do i=1,1000
         read(1,*,end=444) a1,i2,i3,i4,i5
         abad(i)=a1
         if(i4.eq.0.and.i5.eq.0) then
            nbadwa=nbadwa+1
            wbadla(nbadwa)=float(i2)
            wbadua(nbadwa)=float(i3)
         else
            nbada=nbada+1
            ixbsa(nbada)=i2
            ixbea(nbada)=i3
            iybsa(nbada)=i4
            iybea(nbada)=i5
         endif
      enddo
 444  continue
      close(1)

      open(unit=13,file="sub.reg",status="unknown")
      open(unit=101,file="list",status="old")
      ifiball=0
      do iall=1,nmax
         read(101,3001,end=333) file1
         print *,file1

      do i=1,130
         if(file1(i:i+3).eq.'exp0') then
            aef=file1(i:i+4)
            goto 567
         endif
      enddo
 567  continue
      do i=1,130
         if(file1(i:i+4).eq.'multi') then
            amp=file1(i+6:i+20)
            goto 568
         endif
      enddo
 568  continue
      do i=1,130
         if(file1(i:i+2).eq.'201') then
            date=file1(i:i+7)
            goto 569
         endif
      enddo
 569  continue
      do i=1,130
         if(file1(i:i+8).eq.'virus0000') then
            shot=file1(i+9:i+11)
            goto 570
         endif
      enddo
 570  continue

c- check to see if the file exists, and if not write out zeros
         open(unit=1,file=file1,err=111)
         goto 112
 111     continue
         do i=1,nfib
            ifiball=ifiball+1
            do j=1,nw
               xdout(j,ifiball)=0.
            enddo
         enddo
         close(1)
         goto 113
 112     continue
         close(1)
               

      file3="/work/00115/gebhardt/maverick/getampnorm/all/multi_"
     $     //amp//".norm"

c - get the relative frame normalization

      xrelnorm=1.0
      open(unit=1,file='normexp.out',status='old',err=365)
      do i=1,3
         read(1,*,end=366) ae,x2,x3
         if(ae.eq.aef) xrelnorm=x3
      enddo
 366  continue
      close(1)
 365  continue

c- getting the amp2amp

      open(unit=1,file=file3,status='old',err=678)
      na=0
      do i=1,narrm
         read(1,*,end=677) x1,x2
         na=na+1
         wa(na)=x1
         fa(na)=x2
      enddo
 677  continue
      close(1)
      goto 679
 678  continue
      na=3
      wa(1)=3500.
      fa(1)=1.
      wa(2)=4500.
      fa(2)=1.
      wa(3)=5500.
      fa(3)=1.
      write(*,*) "Amp Norm does not exist: ",file3
 679  continue

      call geti(file1,2,xflag,ncol,nrow,ier2)    ! flagged pixels
      call geti(file1,3,x2dsub,ncol,nrow,ier3)   ! 2d sky-sub
      call geti(file1,12,xdw,ncolo,nrow,ier12)   ! wavelength
      call geti(file1,13,xtrace,ncol,nrow,ier13) ! trace
      call geti(file1,14,xdf,ncol,nrow,ier14)    ! f2f
      call geti(file1,16,xdss,ncol,nrow,ier16)   ! sky-subtracted
      call geti(file1,17,xsky,ncol,nrow,ier17)   ! sky
      call geti(file1,18,xderr,ncol,nrow,ier18)  ! error frame

      do ifib=1,nfib
      ifiball=ifiball+1

      n=0
      do i=1,ncolo
         w=xdw(i,ifib)
         if(w.gt.wmin.and.w.lt.wmax) then
            n=n+1
            waved(n)=w
            iwave(n)=i
         endif
      enddo

c - get the bad pixel list
      nbad=0
      nbadw=0
      do i=1,nbada
         a1=abad(i)
         if(a1.eq.amp) then
            nbad=nbad+1
            ixbs(nbad)=ixbsa(i)
            ixbe(nbad)=ixbea(i)
            iybs(nbad)=iybsa(i)
            iybe(nbad)=iybea(i)
         endif
      enddo
      do i=1,nbadwa
         a1=abad(i)
         if(a1.eq.amp) then
            nbadw=nbadw+1
            wbadl(nbadw)=wbadla(i)
            wbadu(nbadw)=wbadua(i)
         endif
      enddo

c - parameters for empirical error estimate
c      nerrsp=51
c      nerrsp2=11
c      nerrf=9
      nerrsp=51
      nerrsp2=3
      nerrf=13
      nerrh=nint(float(nerrsp)/2.)
      nerrh2=nint(float(nerrsp2)/2.)
      nerrfh=nint(float(nerrf)/2.)
      do i=1,n
         xsp(i)=xdss(iwave(i),ifib)
c - get the errors both spectrally and across fibers
c - spectrally:
         iemin=max(1,i-nerrh)
         iemax=min(n,iemin+nerrsp)
         nin=0
         do ie=iemin,iemax
            if(ie.lt.(i-1).or.ie.gt.(i+1)) then
               nin=nin+1
               xin(nin)=xdss(iwave(ie),ifib)
            endif
         enddo
         call biwgt(xin,nin,xb,xerrsp)
         call moment(xin,nin,ave,adev,sdevs,var,skew,curt)
c - across fibers:
         iemin=max(1,i-nerrh2)
         iemax=min(n,iemin+nerrsp2)
         ifmin=max(1,ifib-nerrfh)
         ifmax=min(112,ifmin+nerrf)
         nin=0
         do ie=iemin,iemax
            do ifc=ifmin,ifmax
               if(ifc.ne.ifib) then
                  nin=nin+1
                  xin(nin)=xdss(iwave(ie),ifc)
               endif
            enddo
         enddo
         call biwgt(xin,nin,xb,xerrf) 
         call moment(xin,nin,ave,adev,sdevf,var,skew,curt)
c         xerr(i)=max(xerrsp,xerrf)
         xerr(i)=max(sdevs,sdevf)
      enddo

c - now get the error from the error frame, if it exists

      ierr2=1
      if(ier18.eq.0) then
         do i=1,n
            xval=xderr(iwave(i),ifib)
            if(xval.gt.0.and.xval.lt.1e10) then
               xerr2(i)=xd(iwave(i),ifib)
            else
               xerr2(i)=0.
            endif   
         enddo
      else
         ierr2=0
         do i=1,n
            xerr2(i)=0.
         enddo
c         print *,"Using local error estimate"
      endif   

c - get rid of the horizontal streaks
      nrowe=3
      nrt=nrowe*2+1
      nbin=11
      nhalf=(nbin-1)/2
      diffmax=10.
      xsmin=11.
      do i=1,n
         trace(i)=xtrace(iwave(i),ifib)
         ixtr=nint(trace(i))
         do j=1,nrt
            nr=ixtr-nrowe-1+j
            xspa(j,i)=x2dsub(iwave(i),nr)
         enddo
      enddo      
      do j=1,nrt
         do i=1,n
            is=max(1,i-nhalf)
            ie=min(n,i+nhalf)
            xsum=0.
            do ip=is,ie
               xsum=xsum+xspa(j,ip)
            enddo
            xspab(j,i)=xsum
         enddo
      enddo
      do i=1,n
         ixtr=nint(trace(i))
         nin=0
         do j=1,nrt
            nin=nin+1
            xin(nin)=xspab(j,i)
         enddo
         call biwgt(xin,nin,xbcol,xscol)
         xscol=max(xsmin,xscol)
         do j=1,nrt
            diff=(xspab(j,i)-xbcol)/xscol
            if(diff.gt.diffmax) then
               nr=nint(trace(i)-nrowe-1+j)
               nbad=nbad+1
               ixbs(nbad)=i
               ixbe(nbad)=i
               iybs(nbad)=nr
               iybe(nbad)=nr
            endif
         enddo
      enddo

c - now get the flagged region by wavelength
      if(nbadw.gt.0) then
         do i=1,nbadw
            do j=1,n
               if(waved(j).gt.wbadl(i).and.waved(j).lt.wbadu(i)) then
                  ibad=iwave(j)
                  jbad=nint(xtrace(ibad,ifib))
                  nbad=nbad+1
                  ixbs(nbad)=ibad
                  ixbe(nbad)=ibad
                  iybs(nbad)=jbad
                  iybe(nbad)=jbad
               endif
            enddo
         enddo
      endif

c - now add the extra bad pixel list
      if(nbad.gt.0) then
         do ibad=1,nbad
            do ix=ixbs(ibad),ixbe(ibad)
               do iy=iybs(ibad),iybe(ibad)
                  xflag(ix,iy)=0.
               enddo
            enddo
         enddo
      endif

c      radflag=2.0
      radflagx=3.0
      radflagy=4.0
      do i=1,n
         xftf=xdf(iwave(i),ifib)
         if(xftf.gt.0.1) then
            xsp(i)=xsp(i)/xdf(iwave(i),ifib)
            xerr(i)=xerr(i)/xdf(iwave(i),ifib)
            f2f(i)=xdf(iwave(i),ifib)
            sky(i)=xsky(iwave(i),ifib)

c - flag if anything nearby is below 0
            trace(i)=xtrace(iwave(i),ifib)
            xfpos=float(iwave(i))
            yfpos=trace(i)
            ixmin=max(1,nint(xfpos-radflagx))
            ixmax=min(n,nint(xfpos+radflagx))
            iymin=max(1,nint(yfpos-radflagy))
            iymax=min(n,nint(yfpos+radflagy))
            flag(i)=xflag(iwave(i),nint(trace(i)))
            do ix=ixmin,ixmax
               do iy=iymin,iymax
                  if(xflag(ix,iy).le.0) flag(i)=-1
               enddo
            enddo
         else
            xftf=1.
            xsp(i)=xsp(i)
            xerr(i)=100.
            f2f(i)=xftf
            sky(i)=0.
            flag(i)=-1
         endif
      enddo

c      open(unit=11,file='out.sp',status='unknown')
      do i=1,nw
         call xlinint(wave(i),n,waved,xsp,yp)
         call xlinint(wave(i),n,waved,xerr,yerr)
         call xlinint(wave(i),n,waved,xerr2,yerr2)
         call xlinint(wave(i),n,waved,f2f,yf2f)
         call xlinint(wave(i),n,waved,sky,ysky)
         call xlinint(wave(i),n,waved,trace,ytrace)
         call xlinint(wave(i),n,waved,flag,yflag)
         call xlinint(wave(i),ntp,wtp,tp,ytp)
         call xlinint(wave(i),na,wa,fa,yfp)
c - include the relative dither normalization first:
         if(xrelnorm.gt.0) yp=yp/xrelnorm
         if(xrelnorm.gt.0) yerr=yerr/xrelnorm
c         ytp=ytp*xrelnorm
         yf2f=max(yf2f,0.)
         yfrac=yfp*ytp
c         if(yerr2.gt.5000.or.yerr.gt.5000) yflag=0
         if(yerr2.gt.5000) yflag=0
         if(yerr.gt.1.e5) yflag=0
         if(wave(i).lt.waved(3).or.wave(i).gt.waved(n-3)) yflag=0
         if(yfp.eq.0) yflag=0
c         if(yfp.gt.0.and.yflag.gt.0.and.
         if(yflag.gt.0.and.yp.lt.1e7.and.yp.gt.-1e3) then
c     $        yp.lt.1e5.and.yp.gt.-1e3) then
c - yerr2 is propagation of error
c   ysky is the sky value
c   xrs is 5% of the sky value
c   yerrn is the max of yerr2 or (yerr2+0.05*sky)/1.3
c   or yerrn is the quad sum
c   yerr is the empirical local estimate
            xph=yerr2
            xct=ysky
            xrs=0.04*xct
            yerrn=sqrt(xph*xph+xrs*xrs)
c            yerr2=yerrn
            yerrm=yerr
            if(xph.gt.1) then
               yerr=yerrn
               yerrm=yerrn
            endif
            yerrm=yerrm*convfac/wave(i)/ytp/yfp
c            write(11,1101) wave(i),yp/yfp,yp*convfac/wave(i)/ytp/yfp,
c     $           yfp,ytp,yf2f,yerr2/yfp,yerr/yfp,yerrm
            xdout(i,ifiball)=yp/yfp
         else
c            write(11,1101) wave(i),0.,0.,0.,ytp,yf2f,0.,0.,0.
            xdout(i,ifiball)=0.
         endif
      enddo
c      close(11)
      enddo ! end of fiber loop
 113  continue
      ix=10
      iy=ifiball-56
      write(13,1301) ix,iy,amp
      enddo ! end of file loop
 333  continue
      close(1)
      close(13)
 1301 format("# text("i2,", ",i5,") textangle=90 text={"a14"}")

      naxis=2
      naxes(1)=nw
      naxes(2)=ifiball
      iblock=1
      igc=0
      ier=0

      call ftinit(50,'sub.fits',iblock,ier)
      call ftphps(50,-32,naxis,naxes,ier)
      if(ier.ne.0) then
         print *,'Error in output file ',ier
      endif
      print *,naxes(1),naxes(2)
      call ftp2de(50,igc,1036,naxes(1),naxes(2),xdout,ier)
      call ftclos(50,ier)


 3001 format(a130)

 1101 format(1x,f8.2,1x,f12.3,1x,1pe13.4,3(1x,0pf7.3),2(1x,f9.2),
     $     1x,1pe13.4)
c 1101 format(1x,f8.2,2(1x,f12.2),3(1x,f7.3),2(1x,f9.2),
c     $     1x,f12.2)
 706  continue
      end

      subroutine geti(file1,iext,xd,ncol,nrow,ier)
      parameter (narrm=3000)
      real xd(narrm,narrm)
      integer naxes(2)
      character file1*130
      logical simple,extend,anyf

      im1=0
      ier=0
      iread=0
      im1=50
c      call ftgiou(im1,ier)
      call ftopen(im1,file1,iread,iblock,ier)
      call ftmahd(im1,iext,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      ncol=naxes(1)
      nrow=max(1,naxes(2))
      call ftg2de(im1,igc,0.,narrm,ncol,nrow,xd,anyf,ier)
      call ftclos(im1,ier)
      return
      end

      subroutine xlinint(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            if(y(j).eq.0.or.y(j+1).eq.0) yp=0.
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.gt.x(n)) yp=y(n)
c      if(xp.lt.x(1)) yp=0.
c      if(xp.gt.x(n)) yp=0.
      return
      end

      SUBROUTINE moment(data,n,ave,adev,sdev,var,skew,curt)
      INTEGER n
      REAL adev,ave,curt,sdev,skew,var,data(n)
      INTEGER j
      REAL p,s,ep
c      if(n.le.1)pause 'n must be at least 2 in moment'
      s=0.
      do 11 j=1,n
        s=s+data(j)
11    continue
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      ep=0.
      do 12 j=1,n
        s=data(j)-ave
        ep=ep+s
        adev=adev+abs(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
12    continue
      adev=adev/n
      var=(var-ep**2/n)/(n-1)
      sdev=sqrt(var)
      if(var.ne.0.)then
        skew=skew/(n*sdev**3)
        curt=curt/(n*var**2)-3.
      else
c        pause 'no skew or kurtosis when zero variance in moment'
      endif
      return
      END
