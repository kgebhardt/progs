
      parameter(nmax=10000)
      real w(nmax),wd(nmax),sd(nmax),sumsp(nmax)
      real w1(nmax),w2(nmax),sn(nmax),sumsn(nmax)
      integer ns(nmax)
      character file1*80
      
      ws=3490.
      we=5510.
      wbin=100.
      nbin=0
      do i=1,1000
         wnew=ws+float(i-1)*wbin
         if(wnew.lt.we) then
            nbin=nbin+1
            w(nbin)=wnew+wbin/2.
            w1(nbin)=wnew
            w2(nbin)=wnew+wbin
c            print *,nbin,w(nbin),w1(nbin),w2(nbin)
         else
            goto 766
         endif
      enddo
 766  continue

      do i=1,nbin
         sumsp(i)=0.
         sumsn(i)=0.
         ns(i)=0.
      enddo

      open(unit=1,file='list',status='old')
      do i=1,1000
         read(1,*,end=666) file1
         open(unit=2,file=file1,status='old')
         n=0
         do j=1,nmax
c            read(2,*,end=667) x1,x2
            read(2,*,end=667) x1,x2,x3,x4,x5,x6,x7,x8,x9
            n=n+1
            wd(n)=x1
c            sd(n)=x2
            sd(n)=x4
            sn(n)=x9
         enddo
 667     continue
         close(2)
         do j=1,nbin
            do k=1,n
               if(wd(k).gt.w1(j).and.wd(k).le.w2(j)) then
                  sumsp(j)=sumsp(j)+sd(k)
                  sumsn(j)=sumsn(j)+sn(k)
                  ns(j)=ns(j)+1
               endif
            enddo
         enddo
      enddo
 666  continue
      close(1)

      open(unit=11,file='sumspec.out',status='new')
      do i=1,nbin
         if(ns(i).gt.0) then
            xnorm=sumsn(i)/float(ns(i))
         else
            xnorm=0.
         endif
         write(11,*) w(i),sumsp(i),xnorm
      enddo
      close(11)

      end
