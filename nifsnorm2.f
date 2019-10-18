
      real xin(100000),xw(2),xf(2)

      open(unit=1,file='in',status='old')
      
      w1=22500.
      w2=22900.
      xn1=1.
      xn2=0.75

      xw(1)=22500.
      xw(2)=24000.
      xf(1)=1.03
      xf(2)=0.75

      nb=0
      do i=1,100000
         read(1,*,end=666) x1,x2
         if(x1.gt.w1.and.x1.lt.w2) then
            nb=nb+1
            xin(nb)=x2
         endif
      enddo
 666  continue
      call biwgt(xin,nb,xb,xs)
      rewind(1)

      open(unit=11,file='out',status='unknown')

      do i=1,100000
         read(1,*,end=667) x1,x2
         call xlinint(x1,2,xw,xf,xn)
         write(11,*) x1,x2/xb*xn
      enddo
 667  continue
      close(1)
      close(11)

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
