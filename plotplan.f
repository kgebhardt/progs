
      parameter(nmax=20000)
      real x(nmax),y(nmax),ye(nmax),ya(nmax,100),xin(100)
      real yel(nmax),yeu(nmax),ydiff(nmax),ydiffn(nmax)
      real xl(2),yl(2),xtril(10)
      character file1*80,file2*80,c1*3,ctri(10)*4

      xtri=121./365.
      xtril(1)=40./365.+2018.
      xtril(2)=xtril(1)+xtri
      xtril(3)=xtril(2)+xtri
      xtril(4)=xtril(3)+xtri
      xtril(5)=xtril(4)+xtri
      xtril(6)=xtril(5)+xtri
      ctri(1)="18-1"
      ctri(2)="18-2"
      ctri(3)="18-3"
      ctri(4)="19-1"
      ctri(5)="19-2"
      ctri(6)="19-3"

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

c - diff from Jan 1 to Sept 1 is 122 days
c   let's center on Feb 1, 2018, which is mjd=58150.3
      x0=58130.3
      xt1=58211.3-x0
      xt2=58331.2-x0
      xt3=58452.2-x0
      xt4=xt3+121.
      xt5=xt4+121.
      xt6=xt5+121.
      
      xt1=xt1/365.+2018
      xt2=xt2/365.+2018
      xt3=xt3/365.+2018
      xt4=xt4/365.+2018
      xt5=xt5/365.+2018
      xt6=xt6/365.+2018

      xmin=0.+2018
      xmax=2020.
      ymin=0.
      ymax=2300.
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('Year','Cumulative Number','')

      open(unit=1,file='inall',status='old')
      n=0
      do i=1,100000
         read(1,*,end=667) x1
         x1=x1-7
         if(x1.gt.0) then
            n=n+1
            x(n)=x1
            x(n)=x(n)/365.+2018.
            xcount=float(n)/9.
            y(n)=xcount
         endif
      enddo
 667  continue
      close(1)
      call pgslw(4)
      call pgline(n,x,y)
      open(unit=11,file='out',status='unknown')
      do i=1,n
         write(11,*) x(i),y(i)
      enddo
      close(11)
c      open(unit=1,file='fs2018.plan',status='old')
c      n=0
c      do i=1,nmax
c         read(1,*,end=555) x1,x2
c         n=n+1
c         x(n)=x1
c         y(n)=x2-8
c      enddo
c 555  continue
c      close(1)
c      call pgsci(3)
c      call pgline(n,x,y)
c      call pgsci(1)

      open(unit=1,file='dex.dat',status='old')
      n=0
      do i=1,100000
         read(1,*,end=666) x1,x2
         n=n+1
         x(n)=x1-x0
         if(x1.le.x0) xn0=float(n)
         y(n)=float(n)
      enddo
 666  continue
      close(1)
      do i=1,n
         y(i)=y(i)-xn0
      enddo
      call pgsci(2)
      do i=1,n
         x(i)=x(i)/365.+2018
         y(i)=y(i)+90.
      enddo
      call pgline(n,x,y)
      do i=1050,n
         y(i)=y(i)+90.
      enddo
      call pgsls(2)
c      call pgline(n,x,y)

      call pgsci(1)
      call pgsls(4)
      yl(1)=ymin
      yl(2)=ymax
      xl(1)=xt1
      xl(2)=xl(1)
      call pgline(2,xl,yl)
      xl(1)=xt2
      xl(2)=xl(1)
      call pgline(2,xl,yl)
      xl(1)=xt3
      xl(2)=xl(1)
      call pgline(2,xl,yl)
      xl(1)=xt4
      xl(2)=xl(1)
      call pgline(2,xl,yl)
      xl(1)=xt5
      xl(2)=xl(1)
      call pgline(2,xl,yl)
      xl(1)=xt6
      xl(2)=xl(1)
      call pgline(2,xl,yl)

      call pgsls(1)
      call pgsci(1)
      call pgsch(0.8)
      do i=1,6
         call pgptxt(xtril(i),50.,0.0,0.5,ctri(i))
      enddo

      call pgend

      end
