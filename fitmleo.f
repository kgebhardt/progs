c -- inverts the Jean's equation given vel disp and surface profile

      parameter(nmax=4000,np=200,mm=2,nwk=nmax+6*(nmax*mm+1),mm2=mm*2)

      real f(nmax),rfit(np),xnew(np)
      real fv(nmax),rvfit(np),phi(np)
      real xvnew(np),xm(np),rm(np)
      real fn1(np),fvn1(np),fn(np),yr1(np),yr0(np)
      real dnewdr(np),dv2dr(np),vr1(np),fvn(np),vsql(np)
      real ynewl(np),rl(nmax),rvl(nmax)
      real xnewl(np),xvnewl(np),vsqln(np),dens(np),xmtl(np)
      real*8 y(nmax),c(nmax),wk(nwk),q(mm2),wx(nmax)
      real*8 drv(nmax),cf(nmax),cfv(nmax)
      real*8 dr(nmax),dynewl(nmax),dxm(nmax),drm(nmax),splder,val
      character filein*40,type*3
      external funcf,funcv,func1,func2

      common /ffunc/ dr,cf,m,n
      common /vfunc/ drv,cfv,mv,nd
      common /cfunc/ drm,c,facrho,nall

      parameter (pi=3.141593e0)

      data big,eps,rlim,istepmax /1.e30,1.e-4,1.01,6/
      data gmass /0.8157/

      isigstop=1

c - tidal radii at 432 and 347                                                                                                                                    
      fcon2=3.9
      fcon3=6.8

 1069 format(a)

c -- Let's get the data

 1    call qc1('Input file with fluxes ','fitm.def',filein)
      open(unit=1,file=filein,status='old',err=1)

 2    call qc1('Input file with dispersion ','fitm.def',filein)
      open(unit=2,file=filein,status='old',err=2)

      call qr3('Distance(kpc), radmin, radmax(arcm) ','fitm.def',
     $     distkpc,radmin,radmax)
      amtopc=distkpc*.29
      facm=log10(amtopc*232.8)
      facrho=log10(1./4./pi/amtopc**3)

c      call qi1('Plot Flux stuff ? (1-yes) ','fitm.def',iplot)
      iplot=1
      call savdef

      n=0
      rmin=big
      rmax=-big
      fmin=big
      fmax=-big

      do i=1,nmax
         read(1,*,end=666) x1,x2
         n=n+1
         rl(n)=log10(x1/60.)
c         rl(n)=x1-log10(60.)
         f(n)=log10(10**(-0.4*(x2-4.83))*(20626.5)**2)
c         f(n)=x2
         y(n)=dble(f(n))
         dr(n)=dble(rl(n))
         wx(n)=1.d0
         rmin=min(rmin,rl(n))
         rmax=max(rmax,rl(n))
         fmin=min(fmin,f(n))
         fmax=max(fmax,f(n))
      enddo
 666  continue

      call qc1('Graphics device/type ','fitm.def',type)
c      call pgbegin(0,type,2,2)
      call pgbegin(0,type,2,1)
      call pgscf(2)
      call pgsch(1.9)

      fbit=(fmax-fmin)/10.
      fmin=fmin-fbit
      fmax=fmax+fbit
      rbit=(rmax-rmin)/10.
      rmaxp=rmax+rbit

      if(iplot.eq.1 ) then
         call pgenv(rmin,rmaxp,fmin,fmax,0,30)
         call pgpoint(n,rl,f,21)
         call pglabel('R','\gS','\gS vs. Radius')
      endif

c -- smoothing the flux

      call qd1('Enter val for f ','fitm.def',val)
      call savdef
      md=3
      if(val.eq.0.) md=2
      m=2
      mv=2
      
      call gcvspl(dr,y,nmax,wx,1.d0,m,n,1,md,val,cf,nmax,wk,ier)
      if(ier.ne.0) print *,'ier= ',ier,' in f calc'

c      open(unit=11,file='sigma.out',status='unknown')
c      write(11,*) m,n,nmax,mm2
c      write(11,*) dr,cf,q
c      close(11)
      
      do i=1,np
         rfit(i)=rmin+(rmax-rmin)/float(np-1)*float(i-1)
         jin=i
         fn(i)=sngl(splder(0,m,n,dble(rfit(i)),dr,cf,jin,q))
         fn1(i)=sngl(splder(1,m,n,dble(rfit(i)),dr,cf,jin,q))
c         print *,10**rfit(i),fn(i)
c         write(11,*) 10**rfit(i),10**fn(i)
      enddo
c      close(11)

      if(iplot.eq.1) call pgline(np,rfit,fn)

      nd=0
      rvmin=big
      rvmax=-big
      fvmin=big
      fvmax=-big

c -- getting the flux x v^2

      do i=1,nmax
         read(2,*,end=667) x1,x2
         x1=log10(x1/60.)
c         x2=1
         if(x1.gt.log10(radmin)) then
            nd=nd+1
            do j=1,nd-1
               if(x1.eq.rvl(j)) then
                  x1=rvl(j)+5.e-6
                  nchng=nchng+1
                  print *,'Changed ',x1-5.e-6,x1,nchng
               endif
            enddo
            fp=sngl(splder(0,m,n,dble(x1),dr,cf,n/2,q))
            rvl(nd)=x1
            fv(nd)=2.*log10(x2)+fp
            drv(nd)=dble(rvl(nd))
            y(nd)=dble(fv(nd))
            wx(nd)=1.d0
c            print *,x1,fv(nd)
            rvmin=min(rvmin,rvl(nd))
            rvmax=max(rvmax,rvl(nd))
            fvmin=min(fvmin,fv(nd))
            fvmax=max(fvmax,fv(nd))
         endif
      enddo
 667  continue

      fbit=(fvmax-fvmin)/10.
      fvmin=fvmin-fbit
      fvmax=fvmax+fbit
      rbit=(rvmax-rvmin)/10.
      rvmaxp=rvmax+rbit

      if(iplot.eq.1) then
         call pgenv(rvmin,rvmaxp,fvmin,fvmax,0,30)
         call pgpoint(nd,rvl,fv,21)
         call pglabel('R','\gS x \gs\U2\D\Dp\U',
     $        '\gS x \gs\U2\D\Dp\U vs. Radius')
      endif

      call qd1('Enter val for fv ','fitm.def',val)
      md=3
      if(val.eq.0.) md=2
      call gcvspl(drv,y,nmax,wx,1.d0,mv,nd,1,md,val,cfv,nmax,wk,ier)
      if(ier.ne.0) print *,'ier= ',ier,' in fv calc'

      do i=1,np
         rvfit(i)=rvmin+(rvmax-rvmin)/float(np-1)*float(i-1)
         jin=i
         fvn(i)=sngl(splder(0,mv,nd,dble(rvfit(i)),drv,cfv,jin,q))
         fvn1(i)=sngl(splder(1,mv,nd,dble(rvfit(i)),drv,cfv,jin,q))
      enddo

      if(iplot.eq.1) call pgline(np,rvfit,fvn)


c - tidal radii at 432 and 347
      rvna=432.
      rvnb=347.
      rstarta=rvna
      rstartb=rvnb
      rend=600.
      nstep=1000
      suma=0.
      do i=1,nstep
         rstep=rstarta+float(i-1)*(rend-rstarta)/float(nstep-1)
         fvna=sngl(splder(0,mv,nd,dble(log10(rstep/60.)),drv,cfv,jin,q))
         rstep=rstartb+float(i-1)*(rend-rstartb)/float(nstep-1)
         fvnb=sngl(splder(0,mv,nd,dble(log10(rstep/60.)),drv,cfv,jin,q))
         suma=suma+10**fvna
         sumb=sumb+10**fvnb
      enddo
      suma=suma/float(nstep)
      sumb=sumb/float(nstep)
      print *,suma,sumb

c- getting the deprojected luminosity density

      open(unit=12,file='nu.out',status='unknown')

      xmax=-big
      xmin=big
      ximax=rlim*rl(n)

      factint=2.303/pi/amtopc

      do i=1,np
         s=0.
         do j=1,istepmax
            sold=s
            call midsql2(funcf,rfit(i),rfit(i),ximax,s,j)
            if(abs((s-sold)/s).lt.eps) goto 765
         enddo
 765     continue
         xnew(i)=-s*factint
         xmax=max(xmax,xnew(i))
         xmin=min(xmin,xnew(i))
         write(12,*) 60.*10**rfit(i),xnew(i)
      enddo
      close(12)

c- getting the deprojected density x variance

      xvmax=-big
      xvmin=big
      ximax=rlim*rvl(nd)

      factint=2.303/pi/amtopc
      
      do i=1,np-1
         s=0.
         do j=1,istepmax
            sold=s
            call midsql2(funcv,rvfit(i),rvfit(i),ximax,s,j)
c            print *,s,sold
            if(abs((s-sold)/s).lt.eps) goto 766
         enddo
 766     continue
         xvnew(i)=-s*factint
c         if(xvnew(i).eq.0) xvnew(i)=xvnew(i-1)
c         print *,i,xvnew(i)
         xvmax=max(xvmax,xvnew(i))
         xvmin=min(xvmin,xvnew(i))
c         print *,i,xvnew(i),10**rvfit(i),j
      enddo
      xvnew(np)=xvnew(np-1)

      xmin=big
      xmax=-big
      xvmin=big
      xvmax=-big

      do i=1,np
         xnewl(i)=log10(xnew(i))
         xvnewl(i)=log10(xvnew(i))
         xvmax=max(xvmax,xvnewl(i))
         xvmin=min(xvmin,xvnewl(i))
         xmax=max(xmax,xnewl(i))
         xmin=min(xmin,xnewl(i))
      enddo

      xbit=(xmax-xmin)/10.
      xmin=xmin-xbit
      xmax=xmax+xbit

      xbit=(xvmax-xvmin)/10.
      xvmin=xvmin-xbit
      xvmax=xvmax+xbit

      if(iplot.eq.1) then
         call pgenv(rmin,rmaxp,xmin,xmax,0,30)
         call pgline(np,rfit,xnewl)
         call pglabel('R','\gn','\gn vs. Radius')

         call pgenv(rvmin,rvmaxp,xvmin,xvmax,0,30)
         call pgline(np,rvfit,xvnewl)
         call pglabel('R','\gn x V\U2\D',
     $        '\gn x V\U2\D vs. Radius')
      endif

      if(isigstop.eq.1) goto 555

c- getting deprojected variance

      vsqmin=big
      vsqmax=-big

      nv=0
      do i=1,np
         do j=1,np-1
            if(rvfit(i).ge.rfit(j).and.rvfit(i).lt.rfit(j+1)) then
               xnewp=xnew(j)+((xnew(j+1)-xnew(j))
     $              /(rfit(j+1)-rfit(j)))*
     $              (rvfit(i)-rfit(j))
               goto 600
            endif
         enddo
         goto 601
 600     continue
         nv=nv+1
         vsq=xvnew(i)/xnewp
         ynewl(nv)=log10(xnewp)
         vsql(nv)=log10(vsq)
         rvfit(nv)=rvfit(i)
         dr(nv)=dble(rvfit(nv))
         y(nv)=dble(vsql(nv))
         wx(nv)=1.d0
         dynewl(nv)=dble(ynewl(nv))
c         print *,nv,rvfit(nv),ynewl(nv)
         vsqmin=min(vsqmin,vsql(nv))
         vsqmax=max(vsqmax,vsql(nv))
 601     continue
      enddo

      vbit=(vsqmax-vsqmin)/10.
      vsqmin=vsqmin-vbit
      vsqmax=vsqmax+vbit

      call pgenv(rvmin,rvmaxp,vsqmin,vsqmax,0,30)
      call pgline(nv,rvfit,vsql)
      call pglabel('R','V\U2\D','V\U2\D vs. Radius')

      call qd1('Enter val for vsql ','fitm.def',val)
      md=3
      if(val.eq.0) md=2
      call gcvspl(dr,y,nmax,wx,1.d0,m,nv,1,md,val,c,nmax,wk,ier)

      do i=1,nv
         jin=i
         vsqln(i)=sngl(splder(0,m,nv,dble(rvfit(i)),dr,c,jin,q))
         vr1(i)=sngl(splder(1,m,nv,dble(rvfit(i)),dr,c,jin,q))
      enddo
      call pgline(nv,rvfit,vsqln)

      call qd1('Enter val for ynewl ','fitm.def',val)
      md=3
      if(val.eq.0) md=2
      call gcvspl(dr,dynewl,nmax,wx,1.d0,m,nv,1,md,val,c,nmax,wk,ier)

c      open(unit=11,file='nu.out',status='unknown')
c      write(11,*) m,nv,nmax,mm2
c      write(11,*) dr,c,q
c      close(11)

      do i=1,nv
         jin=i
         yr0(i)=sngl(splder(0,m,nv,dble(rvfit(i)),dr,c,jin,q))
         yr1(i)=sngl(splder(1,m,nv,dble(rvfit(i)),dr,c,jin,q))
      enddo

      xmmin=big
      xmmax=-big
      y1min=big
      y1max=-big
      y2min=big
      y2max=-big
      nall=0

      rmax=log10(radmax)

      do i=1,nv
         rmp=rvfit(i)
         if(rmp.lt.rmax) then
            nall=nall+1
            rm(nall)=rmp
            vsqln(nall)=vsqln(i)
            dnewdr(nall)=yr1(i)
            dv2dr(nall)=vr1(i)
            yr0(nall)=yr0(i)
c            print *,rm(nall),yr1(i),vr1(i)
            y1min=min(y1min,dnewdr(nall))
            y1max=max(y1max,dnewdr(nall))
            y2min=min(y2min,dv2dr(nall))
            y2max=max(y2max,dv2dr(nall))
         endif
      enddo

      xm(1)=1.e-10

      do i=1,nall
         xchk=-(10**(rm(i)+vsqln(i)))*(dnewdr(i)+dv2dr(i))
         if(xchk.gt.0.) then
            xm(i)=log10(xchk)+facm
         else
            xm(i)=xm(max(1,i-1))
         endif
         dxm(i)=dble(xm(i))
         drm(i)=dble(rm(i))
         wx(i)=1.d0
         xmmin=min(xmmin,xm(i))
         xmmax=max(xmmax,xm(i))
c         print *,10**rm(i),10**xm(i)
      enddo

      xmbit=(xmmax-xmmin)/10.
      xmmin=xmmin-xmbit
      xmmax=xmmax+xmbit

      call pgenv(rvmin,rmax,xmmin,xmmax,0,30)
      call pgline(nall,rm,xm)
      call pglabel('R','Mass','Mass vs. Radius')

      call qd1('Enter val for mass ','fitm.def',val)
      md=3
      if(val.eq.0) md=2
      call savdef

      call gcvspl(drm,dxm,nmax,wx,1.d0,m,nall,1,md,val,c,nmax,wk,ier)

      xmmax1=10**xm(nall)
      dmin=big
      dmax=-big
      xmmin=big
      xmmax=-big

      open(unit=11,file='fitm.out',status='unknown')
      open(unit=12,file='numden.dat',status='unknown')

      write(11,*) facm,facrho
      write(12,*) 10**yr0(1),xmmax1,amtopc

c -- calculate phi(r)

      facphi=10**(facm+facrho)
      jin=1
      x1=sngl(splder(0,m,nall,drm(1),drm,c,jin,q))
      x2=sngl(splder(1,m,nall,drm(1),drm,c,jin,q))
      rho0=10**x1*x2*10**facrho/3./facphi

      do i=2,nall-1
         call qromb3(func1,rm(1),rm(i),s1)
         call qromb3(func2,rm(i),rlim*rm(nall),s2)
         phi(i)=-2.303*(s1/10**rm(i)+s2)/facphi-rho0/10**rm(i)
      enddo
      phi(1)=phi(2)
      phi(nall)=phi(nall-1)

c -- calculate density

      dens(1)=1.e-10

      open(unit=13,file='dens.out',status='unknown')
      write(13,*) m,nall,nmax,mm2,facrho
      write(13,*) drm,c,q
      close(13)
      open(unit=13,file='ml.out',status='unknown')

      do i=1,nall
         jin=i
         xm0=sngl(splder(0,m,nall,drm(i),drm,c,jin,q))
         call pgpoint(1,rm(i),xm0,17)
         xm1=sngl(splder(1,m,nall,drm(i),drm,c,jin,q))
         xm2=sngl(splder(2,m,nall,drm(i),drm,c,jin,q))
c         xm0=xm(i)
c         if(i.eq.1) then
c            xm1=0.87
c         else
c            xm1=(xm(i)-xm(i-1))/(rm(i)-rm(i-1))
c         endif
         if(xm1.gt.0.) then
            dens(i)=xm0+log10(xm1)-3.*rm(i)+facrho
            slope=xm2/xm1-2.
c            print *,10**rm(i),slope,(dens(i)-dens(i-1))/(rm(i)-rm(i-1))
         else
            dens(i)=dens(max(1,i-1))
         endif   
         xmtl(i)=10**(dens(i)-yr0(i))
         write(11,*) sngl(drm(i)),dens(i),c(i),ynewl(i),yr0(i),dnewdr(i)
         ar=1./gmass/amtopc*10**(xm0-facm-vsqln(i)-2.*rm(i))
         br=-1./amtopc*dv2dr(i)/10**rm(i)
         vmax2=-2.*phi(i)
         xmlow=gmass*10**vsqln(i)/vmax2
         cr=7./2./phi(i)/amtopc*10**(xm0-facm-2.*rm(i))
         write(12,1121) 10**rm(i),ar,br,dens(i),10**xm0/xmmax1,
     $        xmlow,cr,yr0(i)
         write(13,*) 10**rm(i),xmtl(i)
         dmin=min(dmin,dens(i))
         dmax=max(dmax,dens(i))
         xmmin=min(xmmin,xmtl(i))
         xmmax=max(xmmax,xmtl(i))
      enddo

 1121 format(8(1x,f8.4))

      close(11)
      close(12)
      close(13)

      dbit=(dmax-dmin)/10.
      dmin=dmin-dbit
      dmax=dmax+dbit
      xmbit=(xmmax-xmmin)/10.
      xmmin=xmmin-xmbit
      xmmax=xmmax+xmbit

      call pgenv(rvmin,rmax,dmin,dmax,0,30)
      call pgline(nall,rm,dens)
      call pglabel('R','Density','Density vs. Radius')
      call pgenv(rvmin,rmax,xmmin,xmmax,0,10)
      call pgline(nall,rm,xmtl)
      call pglabel('R','(M/L)\DV\U','M/L vs. Radius')

 555  continue
      call pgend

      end
