
      parameter (narrm=3000,nmax=1000)
      real xdss(narrm,narrm),xderr(narrm,narrm),xda2a(narrm,narrm)
      real wave(narrm),v(narrm*narrm)
      integer naxes(2)
      character file1*130,file2*130,file3*130,amp*14,a1*14,ae*5,aef*5
      character date*8,shot*3,abad(nmax)*14
      logical simple,extend,anyf

c      call pgbegin(0,'?',3,3)
c      call pgpap(0.,1.)
c      call pgscf(2)
c      call pgsch(1.4)
c      call pgslw(2)

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

      open(unit=11,file='out',status='unknown')
      open(unit=12,file='out2',status='unknown')
      open(unit=101,file="list",status="old")
      do iall=1,nmax
         read(101,3001,end=333) file1
c         print *,file1

         do i=1,130
            if(file1(i:i+3).eq.'exp0') then
               aef=file1(i:i+4)
               goto 567
            endif
         enddo
 567     continue
         do i=1,130
            if(file1(i:i+4).eq.'multi') then
               amp=file1(i+6:i+20)
               goto 568
            endif
         enddo
 568     continue
         do i=1,130
            if(file1(i:i+2).eq.'201') then
               date=file1(i:i+7)
               goto 569
            endif
         enddo
 569     continue
         do i=1,130
            if(file1(i:i+8).eq.'virus0000') then
               shot=file1(i+9:i+11)
               goto 570
            endif
         enddo
 570     continue

c- check to see if the file exists, and if not write out zeros
         open(unit=1,file=file1,err=111)
         goto 112
 111     continue
         close(1)
         goto 332
 112     continue
         close(1)

         call geti(file1,16,xdss,ncol,nrow,ier19) ! sky-subtracted
         call geti(file1,17,xda2a,ncol,nrow,ier21) ! a2a
         call geti(file1,18,xderr,ncol,nrow,ier20) ! error frame

         nall=0
         ncheck=0
         ncheck0=8
         do ifib=2,nfib-1
            do iw=2,ncol-1
               do ifc=ifib-1,ifib+1
                  do iwc=iw-1,iw+1
                     if(xda2a(ifc,iwc).le.0) goto 777
                  enddo
               enddo
               st=0.
               se=0.
               do ifc=ifib-1,ifib+1
                  do iwc=iw-1,iw+1
                     st=st+xdss(ifc,iwc)
                     se=se+xderr(ifc,iwc)**2
                  enddo
               enddo
               se=sqrt(se)
               if(se.gt.0) then
                  nall=nall+1
                  v(nall)=st/se
                  if(abs(v(nall)).gt.ncheck0) ncheck=ncheck+1
               endif
 777           continue
            enddo
         enddo
         call biwgt(v,nall,xb,xs)
         write(11,1101) ncheck,nall,xs,file1

         if(ncheck.lt.100.and.xs.lt.1.2) then
            do i=1,nall
               if(abs(v(i)).lt.10.0) then
                  write(12,1201) v(i),ncheck,nall,xs,file1
               endif
            enddo
         endif            

 332     continue
      enddo
 333  continue
      close(1)
      close(11)
      close(12)

 3001 format(a130)
 1101 format(i6,1x,i6,1x,f5.2,3x,a120)
 1201 format(f10.3,1x,i6,1x,i6,1x,f5.2,3x,a120)

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
