      function funcf(x,r)
      parameter(nmax=4000,np=200,mm=2,nwk=nmax+6*(nmax*mm+1),mm2=mm*2)
      real*8 dr(nmax),cf(nmax),splder,q(mm2)
      common /ffunc/ dr,cf,m,n
      xdum1=sngl(splder(0,m,n,dble(x),dr,cf,n/2,q))
      xdum2=sngl(splder(1,m,n,dble(x),dr,cf,n/2,q))
      der=(10**xdum1)*xdum2
      den=10**(2.*x)-10**(2.*r)
      if(den.gt.0.) then
         funcf=der/sqrt(den)
      else
         funcf=0.
      endif
      return
      end

      function funcv(x,r)
      parameter(nmax=4000,np=200,mm=2,nwk=nmax+6*(nmax*mm+1),mm2=mm*2)
      real*8 drv(nmax),cfv(nmax),splder,q(mm2)
      common /vfunc/ drv,cfv,mv,nd
      xdum1=sngl(splder(0,mv,nd,dble(x),drv,cfv,nd/2,q))
      xdum2=sngl(splder(1,mv,nd,dble(x),drv,cfv,nd/2,q))
      der=(10**xdum1)*xdum2
      den=10**(2.*x)-10**(2.*r)
      if(den.gt.0.) then
         funcv=der/sqrt(den)
      else
         funcv=0.
      endif
      return
      end

      function func1(x)
      parameter(nmax=4000)
      real*8 c(nmax),dr(nmax),q(4),splder
      common /cfunc/ dr,c,facrho,n
      x1=sngl(splder(0,2,n,dble(x),dr,c,n/2,q))
      x2=sngl(splder(1,2,n,dble(x),dr,c,n/2,q))
      func1=10**x1*x2*10**facrho
      return
      end

      function func2(x)
      parameter(nmax=4000)
      real*8 c(nmax),dr(nmax),q(4),splder
      common /cfunc/ dr,c,facrho,n
      x1=sngl(splder(0,2,n,dble(x),dr,c,n/2,q))
      x2=sngl(splder(1,2,n,dble(x),dr,c,n/2,q))
      func2=10**x1*x2/10**x*10**facrho
      return
      end

      function func3(x)
      parameter(np=200)
      real rm(np),dens(np)
      common /cfunc3/ rm,dens,nall,xin
      call xlinint(log10(x),nall,rm,dens,d)
      func3=2.*10**d*x/sqrt(x*x-xin*xin)
c      print *,x,xin,d,func3
c      read *
      return
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


      SUBROUTINE midsql2(funk,r,aa,bb,s,n)
      INTEGER n
      REAL aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL ddel,del,sum,tnm,x,func,a,b
      func(x,r)=2.*x*funk(aa+x**2,r)
      b=sqrt(bb-aa)
      a=0.
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b),r)
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x,r)
          x=x+ddel
          sum=sum+func(x,r)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END

      SUBROUTINE qromb3(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-4, JMAX=15, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      write(*,*) '***WARNING: too many steps in qromb3***'
      END
      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) print *,'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

c----------------------------------------------------------------------------- 
      subroutine biwgt (x,n,xbl,xbs)
c----------------------------------------------------------------------------- 

c  This subroutine calculates an estimator of the location and scale
c  of the data set X.  The scale uses the Biweight function in the
c  general formula of "A-estimators." This formula is given on page 416
c  in UREDA (formula 4). The BIWEIGHT function is given by:
c
c                                  u((1-u*u)**2)     abs(u) <= 1
c                         f(u) =
c                                  0                 abs(u) >  1
c
c  where u is defined by
c
c                         u = (X(I) - Med) / (c * Mad)  .
c
c  Med, Mad, and c are the median, the median absolute deviation from
c  the median, and the tuning constant respectively. The tuning
c  constant is a parameter which is chosen depending on the sample
c  size and the specific function being used for the scale estimate.
c  (See page 417 in UREDA).  The tuning constant c is set to 6.0 for 
c  calculation of the location, as recommended by Tukey.
c
c  The biweght location is found iteratively using the formula:
c
c                         T = Med + (sums)
c
c  where (sums) are as given on page 421 in UREDA
c
c  ( some comments on the scale   Here we take c = 9.0.)
c
c  WARNING: This version of the Biweight routine has the side effect of
c  sorting the input array X, and returning the sorted data in place.


      integer n, i
      real x(n), xbl, xbs, xmed, xmad
      real cmad, cmadsq, delta
      real*8 sum1,sum2,t0,t1

      call medmad (x, n, xmed, xmad)

      if (xmad .lt. 1.0e-6 * max(1.0, abs(xmed))) then
         xbl = xmed
         xbs = xmad
         return
      end if

      xbl = xmed
      delta = xmed
      cmad = 6.0 * xmad
      cmadsq = cmad * cmad

      do while (abs(delta) .ge. 1.0e-5 * max(1.0,abs(xbl)))

         sum1 = 0.0
         sum2 = 0.0

         do i=1,n
            t0 = x(i) - xbl
            if (abs(t0) .lt. cmad) then
               t1 = cmadsq - t0 * t0
               t1 = t1 * t1
               sum1 = sum1  + (t0 * t1)
               sum2 = sum2 + t1
            endif
         end do
         
         delta = sum1 / sum2
         xbl = xbl + delta
      end do

      sum1 = 0.0
      sum2 = 0.0
      cmad = 9.0 * xmad
      cmadsq = cmad * cmad
      
      do i=1,n
         t0 = x(i) - xbl
         if (abs(t0) .lt. cmad) then
            t0 = t0 * t0
            t1 = cmadsq - t0
            sum1 = sum1 + (t0 * t1 * t1 * t1 * t1)
            sum2 = sum2 + t1 * (cmadsq - 5.0 * t0)
         endif
      end do

      xbs = n * sqrt(sum1/(n-1.0)) / abs(sum2)

      return
      end



c------------------------------------------------------------------------------
      subroutine medmad (x,n,xmed,xmad)
c------------------------------------------------------------------------------

c  This routine calculates the median and the median absolute deviation of
c  a data set.  The data are supplied in the array X, of length N. The
c  median is returned in XMED, and the median absolute deviation in XMAD.
c
c  If N is odd, the median is the value from the data set that has equal
c  numbers of values greater than it and less than it.  If N is even, the
c  median is the mean of the two central values of the data set.
c
c  The median absolute deviation, as the name implies, is the median of
c  the distribution of absolute values of the deviation of the data about
c  the median.
c
c  WARNING: This routine has the side effect of sorting the input array X.


      integer i, j, k, n, n2
      real x(n), xmed, xmad, xi, xj
      logical even

      even = mod(n,2) .eq. 0
      n2 = n / 2

c  sort the data

      call sort (n,x)

c  calculate the median

      if (even) then
         xmed = 0.5 * (x(n2) + x(n2+1))
      else
         xmed = x(n2+1)
      endif

c  calculate the mad

      i = n2
      j = n2 + 2
      if (even) j = n2 + 1
      xi = xmed - x(i)
      xj = x(j) - xmed

      do k=1,n2
         if (xi .lt. xj) then
            xmad = xi
            if (i .gt. 1) then
               i = i - 1
               xi = xmed - x(i)
            else
               j = j + 1
               xj = x(j) - xmed
            end if
         else
            xmad = xj
            if (j .lt. n) then
               j = j + 1
               xj = x(j) - xmed
            else
               i = i - 1
               xi = xmed - x(i)
            end if
         end if
      end do

      if (even) then
         if (xi .lt. xj) then
            xmad = 0.5 * (xmad + xi)
         else
            xmad = 0.5 * (xmad + xj)
         end if
      end if

      return
      end


c------------------------------------------------------------------------------
      subroutine sort (n,x)
c------------------------------------------------------------------------------
 
c  This routine performs a heapsort of the data in array X, of length N.
c  The data are sorted into ascending order, and returned in the array X.
c
c  Stolen (unabashedly) from NUMERICAL RECIPES.



      integer i, j, k, m, n
      real x(n), temp
 
      m = n / 2 + 1
      k = n

 10   if (m .gt. 1) then
         m = m - 1
         temp = x(m)
      else
         temp = x(k)
         x(k) = x(1)
         k = k - 1
         if (k .eq. 1) then
            x(1) = temp
            return
         endif
      endif

      i = m
      j = m + m

 20   if (j .le. k) then
         if (j .lt. k) then
            if (x(j) .lt. x(j+1)) j = j + 1
         endif
         if (temp .lt. x(j)) then
            x(i) = x(j)
            i = j
            j = j + j
         else
            j = k + 1
         endif
         go to 20
      endif

      x(i) = temp
      go to 10

      end
