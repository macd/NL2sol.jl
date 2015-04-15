      subroutine testr(n, p, x, nfcall, r, uiparm, urparm, ufparm)
c
c     *****parameters.
c
      integer n, p, nfcall, uiparm(1)
      double precision x(p), r(n), urparm(1)
      external ufparm
c
c     ..................................................................
c     ..................................................................
c
c     *****purpose.
c     this routine evaluates  r  for the various test functions in
c        references (1), (2), and (3), as well as for some variations
c        suggested by jorge more (private communication) on some of
c        these test problems (for nex .ge. 30).
c
c     *****parameter description.
c     on input.
c
c        n is the length of r.
c        p is the length of x.
c        x is the point at which the residual vector r is to be
c             computed.
c        nfcall is the invocation count of testr.
c        nex = uiparm(1) is the index of the problem currently being
c             solved.
c        urparm is a user parameter vector (and is ignored).
c        ufparm is a user entry point parameter (and is ignored).
c
c     on output.
c
c        r is the residual vector at x.
c
c     *****application and usage restrictions.
c     these test problems may be used to test least-squares solvers
c     such as nl2sol.  in particular, these problems may be used to
c     check whether  nl2sol  has been successfully transported to
c     a particular machine.
c
c     *****algorithm notes.
c     none
c
c     *****subroutines and functions called.
c     none
c
c     *****references
c     (1) gill, p.e.. & murray, w. (1976), algorithms for the solution
c        of the non-linear least-squares problem, npl report nac71.
c
c     (2) meyer, r.r. (1970), theoretical and computational aspects
c        of nonlinear regression, pp. 465-486 of nonlinear programming,
c        edited by j.b. rosen, o.l.mangasarian, and k. ritter,
c        academic press, new york.
c
c     (3) brown, k.m. (1969), a quadratically convergent newton-
c        like method based upon gaussian elimination,
c        siam j. numer. anal. 6, pp. 560-569.
c
c     *****general.
c
c     this subroutine was written in connection with research
c     supported by the national science foundation under grants
c     mcs-7600324, dcr75-10143, 76-14311dss, and mcs76-11989.
c
c     ..................................................................
c     ..................................................................
c
c  ***  local variables and constants  ***
c
      double precision e1, e2, floatn, ri, r1, r2, t, theta, ti, tim1,
     1             tip1, twopi, t1, t2, u, v, w, z
      double precision ybard(15), ykow(11), ukow(11), yosb1(33),
     1             yosb2(65), ymeyer(16)
      integer i, j, nex, nm1
      double precision expmax, expmin, uftolg
c  ***  intrinsic functions  ***
c/+
      integer mod
      real float
      double precision datan2, dble, dcos, dexp, dlog, dmin1, dsin,
     1                 dsqrt
c/
      external rmdcon
      double precision dfloat, rmdcon
c/6
c      data twopi/6.283185307179586d+0/
c/7
      parameter (twopi=6.283185307179586d+0)
c/
c/6
c/7
      save expmax, expmin, uftolg
c/
      data ybard(1)/1.4d-1/, ybard(2)/1.8d-1/, ybard(3)/2.2d-1/,
     1   ybard(4)/2.5d-1/, ybard(5)/2.9d-1/, ybard(6)/3.2d-1/,
     2   ybard(7)/3.5d-1/, ybard(8)/3.9d-1/, ybard(9)/3.7d-1/,
     3   ybard(10)/5.8d-1/, ybard(11)/7.3d-1/, ybard(12)/9.6d-1/,
     4   ybard(13)/1.34d0/, ybard(14)/2.10d0/, ybard(15)/4.39d0/
      data ykow(1)/1.957d-1/, ykow(2)/1.947d-1/, ykow(3)/1.735d-1/,
     1   ykow(4)/1.600d-1/, ykow(5)/8.44d-2/, ykow(6)/6.27d-2/,
     2   ykow(7)/4.56d-2/, ykow(8)/3.42d-2/, ykow(9)/3.23d-2/,
     3   ykow(10)/2.35d-2/, ykow(11)/2.46d-2/
      data ukow(1)/4.0d0/, ukow(2)/2.0d0/, ukow(3)/1.0d0/,
     1   ukow(4)/5.0d-1/, ukow(5)/2.5d-1/, ukow(6)/1.67d-1/,
     2   ukow(7)/1.25d-1/, ukow(8)/1.0d-1/, ukow(9)/8.33d-2/,
     3   ukow(10)/7.14d-2/, ukow(11)/6.25d-2/
      data yosb1(1)/8.44d-1/, yosb1(2)/9.08d-1/, yosb1(3)/9.32d-1/,
     1   yosb1(4)/9.36d-1/, yosb1(5)/9.25d-1/, yosb1(6)/9.08d-1/,
     2   yosb1(7)/8.81d-1/, yosb1(8)/8.50d-1/, yosb1(9)/8.18d-1/,
     3   yosb1(10)/7.84d-1/, yosb1(11)/7.51d-1/, yosb1(12)/7.18d-1/,
     4   yosb1(13)/6.85d-1/, yosb1(14)/6.58d-1/, yosb1(15)/6.28d-1/,
     5   yosb1(16)/6.03d-1/, yosb1(17)/5.80d-1/, yosb1(18)/5.58d-1/,
     6   yosb1(19)/5.38d-1/, yosb1(20)/5.22d-1/, yosb1(21)/5.06d-1/,
     7   yosb1(22)/4.90d-1/, yosb1(23)/4.78d-1/, yosb1(24)/4.67d-1/,
     8   yosb1(25)/4.57d-1/, yosb1(26)/4.48d-1/, yosb1(27)/4.38d-1/,
     9   yosb1(28)/4.31d-1/, yosb1(29)/4.24d-1/, yosb1(30)/4.20d-1/,
     a   yosb1(31)/4.14d-1/, yosb1(32)/4.11d-1/, yosb1(33)/4.06d-1/
      data yosb2(1)/1.366d0/, yosb2(2)/1.191d0/, yosb2(3)/1.112d0/,
     1   yosb2(4)/1.013d0/, yosb2(5)/9.91d-1/, yosb2(6)/8.85d-1/,
     2   yosb2(7)/8.31d-1/, yosb2(8)/8.47d-1/, yosb2(9)/7.86d-1/,
     3   yosb2(10)/7.25d-1/, yosb2(11)/7.46d-1/, yosb2(12)/6.79d-1/,
     4   yosb2(13)/6.08d-1/, yosb2(14)/6.55d-1/, yosb2(15)/6.16d-1/,
     5   yosb2(16)/6.06d-1/, yosb2(17)/6.02d-1/, yosb2(18)/6.26d-1/,
     6   yosb2(19)/6.51d-1/, yosb2(20)/7.24d-1/, yosb2(21)/6.49d-1/,
     7   yosb2(22)/6.49d-1/, yosb2(23)/6.94d-1/, yosb2(24)/6.44d-1/,
     8   yosb2(25)/6.24d-1/, yosb2(26)/6.61d-1/, yosb2(27)/6.12d-1/,
     9   yosb2(28)/5.58d-1/, yosb2(29)/5.33d-1/, yosb2(30)/4.95d-1/,
     a   yosb2(31)/5.00d-1/, yosb2(32)/4.23d-1/, yosb2(33)/3.95d-1/,
     b   yosb2(34)/3.75d-1/, yosb2(35)/3.72d-1/, yosb2(36)/3.91d-1/,
     c   yosb2(37)/3.96d-1/, yosb2(38)/4.05d-1/, yosb2(39)/4.28d-1/,
     d   yosb2(40)/4.29d-1/, yosb2(41)/5.23d-1/, yosb2(42)/5.62d-1/,
     e   yosb2(43)/6.07d-1/, yosb2(44)/6.53d-1/, yosb2(45)/6.72d-1/,
     f   yosb2(46)/7.08d-1/, yosb2(47)/6.33d-1/, yosb2(48)/6.68d-1/,
     g   yosb2(49)/6.45d-1/, yosb2(50)/6.32d-1/, yosb2(51)/5.91d-1/,
     h   yosb2(52)/5.59d-1/, yosb2(53)/5.97d-1/, yosb2(54)/6.25d-1/,
     i   yosb2(55)/7.39d-1/, yosb2(56)/7.10d-1/, yosb2(57)/7.29d-1/,
     j   yosb2(58)/7.20d-1/, yosb2(59)/6.36d-1/, yosb2(60)/5.81d-1/
      data yosb2(61)/4.28d-1/, yosb2(62)/2.92d-1/, yosb2(63)/1.62d-1/,
     1   yosb2(64)/9.8d-2/, yosb2(65)/5.4d-2/
      data ymeyer(1)/3.478d4/, ymeyer(2)/2.861d4/, ymeyer(3)/2.365d4/,
     1   ymeyer(4)/1.963d4/, ymeyer(5)/1.637d4/, ymeyer(6)/1.372d4/,
     2   ymeyer(7)/1.154d4/, ymeyer(8)/9.744d3/, ymeyer(9)/8.261d3/,
     3   ymeyer(10)/7.030d3/, ymeyer(11)/6.005d3/, ymeyer(12)/5.147d3/,
     4   ymeyer(13)/4.427d3/, ymeyer(14)/3.820d3/, ymeyer(15)/3.307d3/,
     5   ymeyer(16)/2.872d3/
c
      data expmax/0.d0/, uftolg/0.d0/
c
      dfloat(ii) = dble(float(ii))
c
c-----------------------------------------------------------------------
c
      nex = uiparm(1)
      go to (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100,
     1   1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100,
     2   2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 1900, 2100,
     3   2500, 1300, 1400, 1500, 1600), nex
c
c  ***  rosenbrock   ***
 100  r(1) = 1.0d1*(x(2) - x(1)**2)
      r(2) = 1.0d0 - x(1)
      go to 9999
c  ***  helix   ***
 200  theta = datan2(x(2), x(1))/twopi
      if (x(1) .le. 0.d0 .and. x(2) .le. 0.d0) theta = theta + 1.d0
      r(1) = 1.0d1*(x(3) - 1.0d1*theta)
      r(2) = 1.0d1*(dsqrt(x(1)**2 + x(2)**2) - 1.0d0)
      r(3) = x(3)
      go to 9999
c  ***  singular   ***
 300  r(1) = x(1) + 1.0d1*x(2)
      r(2) = dsqrt(5.0d0)*(x(3) - x(4))
      r(3) = (x(2) - 2.0d0*x(3))**2
      r(4) = dsqrt(1.0d1)*(x(1) - x(4))**2
      go to 9999
c  ***  woods   ***
 400  r(1) = 1.0d1*(x(2) - x(1)**2)
      r(2) = 1.0d0 - x(1)
      r(3) = dsqrt(9.0d1)*(x(4) - x(3)**2)
      r(4) = 1.0d0 - x(3)
      r(5) = dsqrt(9.9d0)*(x(2) + x(4) - 2.d0)
      t = dsqrt(2.0d-1)
      r(6) = t*(x(2) - 1.0d0)
      r(7) = t*(x(4) - 1.0d0)
      go to 9999
c  ***  zangwill
 500  r(1) = x(1) - x(2) + x(3)
      r(2) = -x(1) + x(2) + x(3)
      r(3) = x(1) + x(2) - x(3)
      go to 9999
c  ***  engvall   ***
 600  r(1) = x(1)**2 + x(2)**2 + x(3)**2 - 1.0d0
      r(2) = x(1)**2 + x(2)**2 + (x(3) - 2.0d0)**2 - 1.0d0
      r(3) = x(1) + x(2) + x(3) - 1.0d0
      r(4) = x(1) + x(2) - x(3) + 1.0d0
      r(5) = x(1)**3 + 3.0d0*x(2)**2 + (5.0d0*x(3) - x(1) + 1.0d0)**2
     1               - 3.6d1
      go to 9999
c  ***  branin ***
 700  r(1) = 4.0d0*(x(1) + x(2))
      r(2) = r(1) + (x(1) - x(2))*((x(1) - 2.0d0)**2 +
     1       x(2)**2 - 1.0d0)
      go to 9999
c  ***  beale  ***
 800  r(1) = 1.5d0 - x(1)*(1.0d0 - x(2))
      r(2) = 2.25d0 - x(1)*(1.0d0 - x(2)**2)
      r(3) = 2.625d0 - x(1)*(1.0d0 -  x(2)**3)
      go to 9999
c  ***  cragg and levy  ***
 900  r(1) = (dexp(x(1)) - x(2))**2
      r(2) = 1.0d1*(x(2) - x(3))**3
      r(3) = ( dsin(x(3) - x(4)) / dcos(x(3) - x(4)) )**2
      r(4) = x(1)**4
      r(5) = x(4) - 1.0d0
      go to 9999
c  ***  box  ***
 1000 if (expmax .gt. 0.d0) go to 1001
         expmax = 1.999d0 * dlog(rmdcon(5))
         expmin = 1.999d0 * dlog(rmdcon(2))
 1001 if (-expmax .ge. dmin1(x(1), x(2), x(3))) go to 1003
      do 1002 i = 1,10
         ti = -0.1d0*dfloat(i)
         t1 = ti*x(1)
         e1 = 0.d0
         if (t1 .gt. expmin) e1 = dexp(t1)
         t2 = ti*x(2)
         e2 = 0.d0
         if (t2 .gt. expmin) e2 = dexp(t2)
         r(i) = (e1 - e2) - x(3)*(dexp(ti) - dexp(1.0d1*ti))
 1002 continue
      go to 9999
 1003 nfcall = -1
      go to 9999
c  ***  davidon 1  ***
 1100 nm1 = n - 1
      do 1102 i = 1, nm1
         r1 = 0.0d0
         ti = dfloat(i)
         t = 1.d0
         do 1101 j = 1,p
              r1 = r1 + t*x(j)
              t = t*ti
 1101         continue
         r(i) = r1
 1102    continue
      r(n) = x(1) - 1.0d0
      go to 9999
c  ***  freudenstein and roth  ***
 1200 r(1) = -1.3d1 + x(1) - 2.0d0*x(2) + 5.0d0*x(2)**2 - x(2)**3
      r(2) = -2.9d1 + x(1) - 1.4d1*x(2) + x(2)**2 + x(2)**3
      go to 9999
c  ***  watson  ***
 1300  continue
 1400  continue
 1500  continue
 1600 do 1602 i = 1, 29
         ti = dfloat(i)/2.9d1
         r1 = 0.0d0
         r2 = x(1)
         t = 1.0d0
         do 1601 j = 2, p
              r1 = r1 + dfloat(j-1)*t*x(j)
              t = t*ti
              r2 = r2 + t*x(j)
 1601         continue
         r(i) = r1 - r2*r2 - 1.0d0
         if (nex .ge. 33 .and. nex .le. 36) r(i) = r(i) + 10.d0
 1602    continue
      r(30) = x(1)
      r(31) = x(2) - x(1)**2 - 1.0d0
      if (nex .lt. 33 .or. nex .gt. 36) go to 9999
      r(30) = r(30) + 10.d0
      r(31) = r(31) + 10.d0
      go to 9999
c  ***  chebyquad  ***
 1700 do 1701 i = 1,n
 1701    r(i) = 0.0d0
      do 1702 j = 1,n
         tim1 = 1.0d0
         ti = 2.0d0*x(j) - 1.0d0
         z = ti + ti
         do 1702 i = 1,n
              r(i) = r(i) + ti
              tip1 = z*ti -tim1
              tim1 = ti
              ti = tip1
 1702         continue
      floatn = dfloat(n)
      do 1703 i = 1,n
         ti = 0.0d0
         if (mod(i,2) .eq. 0) ti = -1.0d0/dfloat(i*i - 1)
         r(i) = ti - r(i)/floatn
 1703    continue
      go to 9999
c  ***  brown and dennis  ***
 1800  do 1801 i = 1, n
         ti = 0.2d0*dfloat(i)
         r(i) = (x(1) + x(2)*ti - dexp(ti))**2 +
     1             (x(3) + x(4)*dsin(ti) - dcos(ti))**2
 1801    continue
      go to 9999
c  ***  bard  ***
 1900 do 1901 i = 1, 15
         u = dfloat(i)
         v = 1.6d1 - u
         w = dmin1(u,v)
         r(i) = ybard(i) - (x(1) + u/(x(2)*v + x(3)*w))
         if (nex .eq. 30) r(i) = r(i) + 10.d0
 1901    continue
      go to 9999
c  ***  jennrich and sampson  ***
 2000 do 2001 i = 1, 10
         ti = dfloat(i)
         r(i) = 2.0d0 + 2.0d0*ti - (dexp(ti*x(1)) +
     1          dexp(ti*x(2)))
 2001    continue
      go to 9999
c  ***  kowalik and osborne  ***
 2100 do 2101 i = 1, 11
         r(i) = ykow(i) - x(1)*(ukow(i)**2 + x(2)*ukow(i))/(ukow(i)**2 +
     1          x(3)*ukow(i) + x(4))
         if (nex .eq. 31) r(i) = r(i) + 10.d0
 2101    continue
      go to 9999
c  ***  osborne 1  ***
 2200 do 2201 i = 1, 33
         ti = 1.0d1*dfloat(1-i)
         r(i) = yosb1(i) - (x(1) + x(2)*dexp(x(4)*ti) +
     1          x(3)*dexp(x(5)*ti))
 2201    continue
      go to 9999
c  ***  osborne 2  ***
c     ***  uftolg is a machine-dependent constant.  it is just slightly
c     ***  larger than the log of the smallest positive machine number.
 2300 if (uftolg .eq. 0.d0) uftolg = 1.999d0 * dlog(rmdcon(2))
      do 2302 i = 1, 65
         ti = 0.1d0*dfloat(1-i)
         ri = x(1)*dexp(x(5)*ti)
         do 2301 j = 2, 4
              t = 0.d0
              theta = -x(j+4) * (ti + x(j+7))**2
              if (theta .gt. uftolg) t = dexp(theta)
              ri = ri + x(j)*t
 2301         continue
         r(i) = yosb2(i) - ri
 2302 continue
      go to 9999
c  ***  madsen  ***
 2400 r(1) = x(1)**2 + x(2)**2 + x(1)*x(2)
      r(2) = dsin(x(1))
      r(3) = dcos(x(2))
      go to 9999
c  ***  meyer  ***
 2500 do 2501 i = 1, 16
         ti = dfloat(5*i + 45)
         r(i)=x(1)*dexp(x(2)/(ti + x(3))) - ymeyer(i)
         if (nex .eq. 32) r(i) = r(i) + 10.d0
 2501    continue
      go to 9999
c  ***  brown  ***
 2600 continue
 2700 continue
 2800 continue
 2900 t = x(1) - dfloat(n + 1)
      do 2901 i = 2, n
 2901    t = t + x(i)
      nm1 = n - 1
      do 2902 i = 1, nm1
 2902    r(i) = t + x(i)
      t = x(1)
      do 2903 i = 2, n
 2903    t = t * x(i)
      r(n) = t - 1.0d0
      go to 9999
c
 9999 return
c     ..... last card of testr .........................................
      end
