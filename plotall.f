
      parameter(nmax=100000,np=7)
      real x(nmax,np),xmina(np),xmaxa(np)
      integer ip(nmax),iclass(nmax)
      character a1*70,lab(np)*20
      common/cosmo/ h0,xmat,xlam,xr

      open(unit=1,file='in',status='old')

      do i=1,nmax
         do j=1,np
            x(i,j)=0.
         enddo
      enddo

      do i=1,nmax
         read(1,*,end=666) i1,(x(i,j),j=1,np)
         iclass(i)=i1
         n=i
      enddo
 666  continue
      close(1)

      open(unit=1,file='lab.in',status='old')
      do i=1,np
         read(1,*,end=668) a1,x1,x2
         lab(i)=a1
         xmina(i)=x1
         xmaxa(i)=x2
      enddo
 668  continue

      call pgbegin(0,'?',np,np)
      call pgscf(2)
      call pgslw(1)
      call pgpap(0.,1.)

      do i=1,np
         do j=1,np
            if(i.eq.j) then
               call pgpage
               call pgvport(0.,1.,0.,1.)
               call getlim1(i,x,n,xmin,xmax)
c               call pgwindow(xmin,xmax,xmin,xmax)
               call pgwindow(xmina(i),xmaxa(i),xmina(i),xmaxa(i))
               call pgsch(7.0)
c               xcen=(xmax+xmin)/2.
               xcen=(xmaxa(i)+xmina(i))/2.
               call pgptxt(xcen,xcen,0.,0.5,lab(i))
               call pgsch(4.0)
               if(i.eq.np) then
                  call pgvport(0.05,0.95,0.15,0.85)
                  call pgbox('m',0.,0,'',0.,0)
                  call pgvport(0.15,0.85,0.05,0.95)
                  call pgbox('',0.,0,'n',0.,0)
               else
                  call pgvport(0.05,0.95,0.15,0.85)
                  call pgbox('n',0.,0,'',0.,0)
                  call pgvport(0.15,0.85,0.05,0.95)
                  call pgbox('',0.,0,'m',0.,0)
               endif
               call pgsch(1.0)
            else
               call plot(j,i,x,n,iclass,xmina,xmaxa)
            endif
         enddo
      enddo

      call pgend

      end

      subroutine plot(ip,jp,xa,n,iclass,xmina,xmaxa)
      parameter(nmax=100000,np=7)
      real xa(nmax,np),x(nmax),y(nmax)
      real xmina(np),xmaxa(np)
      integer iclass(nmax)
      do i=1,n
         x(i)=xa(i,ip)
      enddo
      do i=1,n
         y(i)=xa(i,jp)
      enddo
      call getlim(n,x,y,xmin,xmax,ymin,ymax)
      call pgpage
      call pgvport (0.05,0.95,0.05,0.95)
c      call pgwindow(xmin,xmax,ymin,ymax)
      call pgwindow(xmina(ip),xmaxa(ip),xmina(jp),xmaxa(jp))
      call pgsch(1.5)
      call pgbox('bct',0.,0,'bct',0.,0)
      call pgsch(7.0)
      itype=17
      do i=1,n
         call pgsch(1.0)
         if(iclass(i).ge.4) then
            call pgsci(4)
            call pgpt(1,x(i),y(i),itype)
         elseif(iclass(i).eq.2) then
            call pgsci(3)
            call pgpt(1,x(i),y(i),itype)
         else
            call pgsch(2.0)
            call pgsci(2)
         endif
      enddo
      call pgsci(1)
      call pgsch(1.0)
      return
      end

      subroutine getlim(n,x,y,xmin,xmax,ymin,ymax)
      real x(n),y(n)
      data big/1.e20/
      xmin=big
      xmax=-big
      ymin=big
      ymax=-big
      do i=1,n
         if(nint(x(i)).ne.9999) then
            xmin=min(xmin,x(i))
            xmax=max(xmax,x(i))
         endif
         if(nint(y(i)).ne.9999) then
            ymin=min(ymin,y(i))
            ymax=max(ymax,y(i))
         endif
      enddo
      xbit=(xmax-xmin)/10.
      xmin=xmin-xbit
      xmax=xmax+xbit
      ybit=(ymax-ymin)/10.
      ymin=ymin-ybit
      ymax=ymax+ybit
      return
      end

      subroutine getlim1(ip,x,n,xmin,xmax)
      parameter(nmax=100000,np=7)
      real x(nmax,np)
      data big/1.e20/
      xmin=big
      xmax=-big
      ymin=big
      ymax=-big
      do i=1,n
         val=x(i,ip)
         if(nint(val).ne.9999) then
            xmin=min(xmin,val)
            xmax=max(xmax,val)
         endif
      enddo
      xbit=(xmax-xmin)/10.
      xmin=xmin-xbit
      xmax=xmax+xbit
      return
      end
