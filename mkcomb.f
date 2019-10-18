
      parameter(nmax=10000)
      real x(nmax),y(nmax),xina(nmax,nmax),xin(nmax),yin(nmax)
      real ye(nmax)
      integer id(nmax),idn(nmax)
      character c1*11,c2*11

      c1="tmp???.dat"
      c2="tmp???.out"
      
      open(unit=1,file='rclist',status='old')

      n=0
      do i=1,nmax
         read(1,*,end=666) i1,x2,x3
         n=n+1
         id(n)=i1
         x(n)=x2
         y(n)=x3
      enddo
 666  continue

      rad=2.
      do i=1,n
         nin=0
         do j=1,n
            dist=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2)
            if(dist.le.rad) then
               nin=nin+1
               idn(nin)=id(j)
            endif
         enddo
         do j=1,nin
            write(c1(4:6),1001) idn(j)
            open(unit=1,file=c1,status='old')
            nw=0
            do k=1,nmax
               read(1,*,end=667) x1,x2,x3,x4,x5,x6,x7,x8
               nw=nw+1
               xina(j,nw)=x2
               ye(k)=x8
            enddo
 667        continue
            close(1)
         enddo
         do k=1,nw
            do j=1,nin
               xin(j)=xina(j,k)
            enddo
            call biwgt(xin,nin,xb,xs)
            yin(k)=xb
         enddo
         write(c2(4:6),1001) id(i)
         open(unit=11,file=c2,status='unknown')
         do k=1,nw
            wave=3470+float(k-1)*2
            write(11,*) wave,yin(k),ye(k)/sqrt(float(nin))
         enddo
         close(11)
      enddo

 1001 format(i3)
      end

      
