
      parameter(nmax=1000)
      real r(nmax),v(nmax),ra(nmax,10),va(nmax,10)
      real xl(2),yl(2),rrh(nmax),rath(nmax),suma(nmax)
      real sumall(nmax,nmax),suml(nmax),sumh(nmax)
      integer iang(10)

      character file1*50,gname*50,file2*50

      na=5
      open(unit=1,file='galint.dat',status='old')
      read(1,*,end=667) gname,file1,rmin,rmax,vmin,vmax,rhalf,icore
c         read(1,*) file2
      close(1)
      open(unit=1,file='intlist',status='old')

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgslw(2)
      call pgsch(1.4)

      read(1,*)
      do iall=1,1000
         read(1,*,end=667) file2
         open(unit=2,file=file2,status='old')
         read(2,*)
         n=0
         ia=1
         diff=1.e10
         do i=1,nmax
            read(2,*,end=666) x1,x2,x3,x4,x5,x6,x7
            if(i.eq.1) x2old=x2
            if(x2.ne.x2old) then
               ia=ia+1
               n=0
            endif
            n=n+1
            ra(n,ia)=log10(x1)
            vr=x3
            vt=x4
            vp=x6
            vprot=x7
            vp=sqrt(vp*vp+vprot*vprot)
            vtang=sqrt((vt*vt+vp*vp)/2.)
c            va(n,ia)=vr/vt
            va(n,ia)=vr/vtang
            x2old=x2
c            print *,n,ia,ra(n,ia),va(n,ia)
         enddo
 666     continue
         close(2)
         nt=n
         iamax=ia
         rm=1e10
         if(iall.eq.1) then
            call pgpage
            call pgvport(0.20,0.95,0.15,0.90)
            call pgwindow(log10(rmin),log10(rmax),vmin,vmax)
            call pgbox('bclnst',0.,0,'bcnst',0.,0)
            call pgsch(1.4)
            call pgmtxt('B',2.5,0.5,0.5,
     $           'Radius (arcsecond)')
            call pglabel('','\gs\Dr\U/\gs\Dt','')
c            call pgmtxt('T',-1.7,0.5,0.5,gname)
            xl(1)=log10(rmin)
            xl(2)=log10(rmax)
            yl(1)=1.
            yl(2)=1.
            call pgline(2,xl,yl)
         endif
         call pgsch(3.0)
         do ia=1,na
            do j=1,nt
               r(j)=ra(j,ia)
               v(j)=va(j,ia)
            enddo
         enddo
         do j=1,nt
            sum=0.
            isum=0
            do ia=1,na
               sum=sum+va(j,ia)
               isum=isum+1
            enddo
            suma(j)=sum/float(isum)
            sumall(iall,j)=suma(j)
            nall=iall
         enddo
      enddo
 667  continue
      do j=1,nt
         avg=0.
         rmax=-1e10
         rmin=1e10
         do i=1,nall
            avg=avg+sumall(i,j)
            rmax=max(rmax,sumall(i,j))
            rmin=min(rmin,sumall(i,j))
         enddo
         avg=avg/float(nall)
         suma(j)=avg
         suml(j)=rmin
         sumh(j)=rmax
         print *,10**r(j),suma(j),suml(j),sumh(j)
      enddo
      call pgsci(1)
      call pgslw(4)
      call pgsls(1)
      call pgline(nt,r,suma)
      call pgsls(4)
      call pgline(nt,r,suml)
      call pgline(nt,r,sumh)

      rewind(1)
      do iall=1,1
         read(1,*,end=767) file2
         open(unit=2,file=file2,status='old')
         read(2,*)
         n=0
         ia=1
         diff=1.e10
         do i=1,nmax
            read(2,*,end=766) x1,x2,x3,x4,x5,x6,x7
            if(i.eq.1) x2old=x2
            if(x2.ne.x2old) then
               ia=ia+1
               n=0
            endif
            n=n+1
            ra(n,ia)=log10(x1)
            vr=x3
            vt=x4
            vp=x6
            vprot=x7
            vp=sqrt(vp*vp+vprot*vprot)
            vtang=sqrt((vt*vt+vp*vp)/2.)
c            va(n,ia)=vr/vt
            va(n,ia)=vr/vtang
            x2old=x2
c            print *,n,ia,ra(n,ia),va(n,ia)
         enddo
 766     continue
         close(2)
         nt=n
         iamax=ia
         rm=1e10
         do ia=1,na
            do j=1,nt
               r(j)=ra(j,ia)
               v(j)=va(j,ia)
            enddo
         enddo
         do j=1,nt
            sum=0.
            isum=0
            do ia=1,na
               sum=sum+va(j,ia)
               isum=isum+1
            enddo
            suma(j)=sum/float(isum)
            sumall(iall,j)=suma(j)
            nall=iall
         enddo
      enddo
 767  continue
      close(1)
      do j=1,nt
         avg=0.
         rmax=-1e10
         rmin=1e10
         do i=1,nall
            avg=avg+sumall(i,j)
            rmax=max(rmax,sumall(i,j))
            rmin=min(rmin,sumall(i,j))
         enddo
         avg=avg/float(nall)
         suma(j)=avg
         suml(j)=rmin
         sumh(j)=rmax
      enddo
      call pgsci(3)
      call pgslw(4)
      call pgsls(2)
c      call pgline(nt,r,suma)

      call pgend
      end

      subroutine getlim(file2,iang,rm)
      parameter(nb=10)
      integer iang(nb)
      character file2*50

      do i=1,nb
         iang(i)=0
      enddo
      rm=0.
      open(unit=12,file=file2,status='old')
      do i=1,1000
         read(12,*,end=666) i1,x2
         iang(i1)=1
         rm=max(rm,x2)
      enddo
 666  continue
      close(12)
      return
      end
