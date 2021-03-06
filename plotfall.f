
      parameter (radtodeg=57.29578)
      real xp(100000),yp(100000),xo(1000)
      character a1*1

      xmin=6.3
      xmax=33.3
      ymin=-0.11
      ymax=0.11
      rat=(ymax-ymin)/(xmax-xmin)
      print *,1./rat
      
      np=10
      xmin0=xmin
      xrange=(xmax-xmin)/float(np)
      xmax0=xmin+xrange
      call pgbegin(0,'?',1,10)
c      call pgpap(0.,1)
      call pgscf(2)
      call pgsch(2.5)
      do ip=1,10
         xmin=xmin0+float(ip-1)*xrange
         xmax=xmin+xrange
         print *,xmin,xmax
c      call pgpap(0.,rat)
         call pgenv(xmin,xmax,ymin,ymax,0,0)
         call pglabel('RA(deg)','DEC','')
         
         open(unit=1,file='j3',status='old')
         do i=1,100000
            read(1,*,end=669) x1,x2,x3,x4,x5,x6,x7,x8
            xp(1)=x1
            yp(1)=x2
            xp(2)=x3
            yp(2)=x4
            xp(3)=x5
            yp(3)=x6
            xp(4)=x7
            yp(4)=x8
            call pgpoly(4,xp,yp)
         enddo
 669     continue
         close(1)
      enddo

      call pgend
      end

