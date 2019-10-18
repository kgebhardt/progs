
      parameter(nmax=10000)
      real w(nmax),wo(nmax,nmax),xin(nmax),yin(nmax)

      character file1*100
      open(unit=1,file='inlist',status='old')

      ntot=0
      do ia=1,nmax
         read(1,*,end=666) file1
         ntot=ntot+1
         open(unit=2,file=file1,status='old')
         n=0
         do i=1,nmax
            read(2,*,end=667) x1,x2
            n=n+1
            w(n)=x1
            wo(ia,n)=x2
         enddo
 667     continue
         close(2)
      enddo
 666  continue
      close(1)

      open(unit=11,file='ata.out',status='new')
      do i=1,n
         do j=1,ntot
            xin(j)=wo(j,i)
         enddo
         call biwgt(xin,ntot,xb,xs)
         write(11,*) w(i),xb,xs,ntot
      enddo
      close(11)

      end
