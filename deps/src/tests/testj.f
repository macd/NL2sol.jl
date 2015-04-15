      subroutine testj(n, p, x, nfcall, j, uiparm, urparm, ufparm)
c
c  ***  parameters  ***
c
      integer n, p, nfcall, uiparm(1)
      double precision x(p), j(n,p), urparm(1)
      external ufparm
c
c     ..................................................................
c     ..................................................................
c
c     *****purpose.
c     this routine evaluates the jacobian matrix  j  for the various
c     test problems listed in references (1), (2), and (3).
c
c     *****parameter description.
c     on input.
c
c        nn is the row dimension of  j  as declared in the calling
c             program.
c        n is the actual number of rows in  j  and is the length of  r.
c        p is the number of parameters being estimated and hence is
c             the length of x.
c        x is the vector of parameters at which the jacobian matrix  j
c             is to be computed.
c        nfcall is the invocation count of  testr  at the time when  r
c             was evaluated at  x.  testr ignores nfcall.
c        r is the residual vector at  x  (and is ignored).
c        nex = uiparm(1) is the index of the problem currently being
c             solved.
c        urparm is a user parameter vector (and is ignored).
c        ufparm is a user entry point parameter (and is ignored).
c        testr is the subroutine that computes  r  (and is ignored).
c
c     on output.
c
c        j is the jacobian matrix at x.
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
c     (1) gill, p.e.; & murray, w. (1976), algorithms for the solution
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
      double precision e, expmin, r2, t, theta, ti, tim1, tip1, tpi,
     1   tpim1, tpip1, twopi, u, uftolg, ukow(11), v, w, z, zero
      integer i, k, nex, nm1
c  ***  intrinsic functions  ***
c/+
      real float
      double precision dble, dcos, dexp, dlog, dmin1, dsin, dsqrt
c/
      external rmdcon
      double precision dfloat, rmdcon
c
c/6
c/6                                                                    t
c      data twopi/6.283185307179586d+0/, zero/0.d+0/
c/7
      parameter (twopi=6.283185307179586d+0, zero=0.d+0)
c/
c/6
c/7
      save expmin, uftolg
c/
      data ukow(1)/4.0d0/, ukow(2)/2.0d0/, ukow(3)/1.0d0/,
     1   ukow(4)/5.0d-1/, ukow(5)/2.5d-1/, ukow(6)/1.67d-1/,
     2   ukow(7)/1.25d-1/, ukow(8)/1.0d-1/, ukow(9)/8.33d-2/,
     3   ukow(10)/7.14d-2/, ukow(11)/6.25d-2/
c  ***  machine dependent constant  ***
      data expmin/0.0d0/, uftolg/0.d0/
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
c  ***  rosenbrock  ***
 100  j(1,1) = -2.0d1*x(1)
      j(1,2) = 1.0d1
      j(2,1) = -1.0d0
      j(2,2) = 0.0d0
      go to 9999
c  ***  helix  ***
 200  t = x(1)**2 + x(2)**2
      ti = 1.d2/(twopi*t)
      j(1,1) = ti*x(2)
      t = 1.d1/dsqrt(t)
      j(2,1) = x(1)*t
      j(3,1) = 0.d0
      j(1,2) = -ti*x(1)
      j(2,2) = x(2)*t
      j(3,2) = 0.d0
      j(1,3) = 1.d1
      j(2,3) = 0.d0
      j(3,3) = 1.d0
      go to 9999
c  ***  singular  ***
 300  do 301 k = 1,4
         do 301 i = 1,4
 301          j(i,k) = 0.d0
      j(1,1) = 1.d0
      j(1,2) = 1.d1
      j(2,3) = dsqrt(5.d0)
      j(2,4) = -j(2,3)
      j(3,2) = 2.d0*(x(2) - 2.d0*x(3))
      j(3,3) = -2.d0*j(3,2)
      j(4,1) = dsqrt(4.d1)*(x(1) - x(4))
      j(4,4) = -j(4,1)
      go to 9999
c  ***  woods  ***
 400  do 401 k = 1,4
         do 401 i = 1,7
 401            j(i,k) = 0.d0
      j(1,1) = -2.d1*x(1)
      j(1,2) = 1.d1
      j(2,1) = -1.d0
      j(3,4) = dsqrt(9.d1)
      j(3,3) = -2.d0*x(3)*j(3,4)
      j(4,3) = -1.d0
      j(5,2) = dsqrt(9.9d0)
      j(5,4) = j(5,2)
      j(6,2) = dsqrt(0.2d0)
      j(7,4) = j(6,2)
      go to 9999
c  ***  zangwill  ***
 500  do 501 k = 1,3
         do 501 i = 1,3
 501            j(i,k) = 1.d0
      j(1,2) = -1.d0
      j(2,1) = -1.d0
      j(3,3) = -1.d0
      go to 9999
c  ***  engvall  ***
 600  j(1,1) = 2.d0*x(1)
      j(1,2) = 2.d0*x(2)
      j(1,3) = 2.d0*x(3)
      j(2,1) = j(1,1)
      j(2,2) = j(1,2)
      j(2,3) = 2.d0*(x(3) - 2.d0)
      j(3,1) = 1.d0
      j(3,2) = 1.d0
      j(3,3) = 1.d0
      j(4,1) = 1.d0
      j(4,2) = 1.d0
      j(4,3) = -1.d0
      t = 2.d0*(5.d0*x(3) - x(1) + 1.d0)
      j(5,1) = 3.d0*x(1)**2 - t
      j(5,2) = 6.d0*x(2)
      j(5,3) = 5.d0*t
      go to 9999
c  ***  branin  ***
 700  j(1,1) = 4.d0
      j(1,2) = 4.d0
      j(2,1) = 3.d0 + (x(1) - 2.d0)*(3.d0*x(1) - 2.d0*x(2) - 2.d0) +
     1   x(2)*x(2)
      j(2,2) = 1.d0 + 2.d0*(2.d0*x(1) - x(2)*x(2)) - (x(1) - x(2))**2
      go to 9999
c  ***  beale  ***
 800  j(1,1) = x(2) - 1.d0
      j(1,2) = x(1)
      j(2,1) = x(2)**2 - 1.d0
      j(2,2) = 2.d0*x(1)*x(2)
      j(3,1) = x(2)**3 - 1.d0
      j(3,2) = 3.d0*x(1)*(x(2)**2)
      go to 9999
c  ***  cragg & levy  ***
 900  do 901 i = 1,5
         do 901 k = 1,4
 901          j(i,k) = 0.d0
      t = dexp(x(1))
      j(1,2) = -2.d0*(t - x(2))
      j(1,1) = -t * j(1,2)
      j(2,2) = 3.0d1*(x(2) - x(3))**2
      j(2,3) = -j(2,2)
      j(3,3) = 2.d0*dsin(x(3) - x(4))/(dcos(x(3) - x(4)))**3
      j(3,4) = -j(3,3)
      j(4,1) = 4.d0*x(1)**3
      j(5,4) = 1.d0
      go to 9999
c  ***  box  ***
 1000 if (expmin .eq. zero) expmin = 1.999d0*dlog(rmdcon(2))
      do 1001 i = 1,10
         ti = -0.1d0*dfloat(i)
         e = zero
         t = x(1)*ti
         if (t .ge. expmin) e = dexp(t)
         j(i,1) = ti*e
         e = zero
         t = x(2)*ti
         if (t .ge. expmin) e = dexp(t)
         j(i,2) = -ti*e
         j(i,3) = dexp(1.d1*ti) - dexp(ti)
 1001    continue
      go to 9999
c  ***  davidon 1  ***
 1100 nm1 = n-1
      do 1101 i = 1,nm1
         ti = dfloat(i)
         t = 1.d0
         do 1101 k = 1,p
              j(i,k) = t
              t = t*ti
 1101         continue
      j(n,1) = 1.d0
      do 1102 k = 2,p
 1102    j(n,k) = 0.d0
      go to 9999
c  ***  freudenstein & roth  ***
 1200 j(1,1) = 1.d0
      j(1,2) = -2.d0 + x(2)*(1.d1 - 3.d0*x(2))
      j(2,1) = 1.d0
      j(2,2) = -1.4d1 + x(2)*(2.d0 + 3.d0*x(2))
      go to 9999
c  ***  watson  ***
 1300 continue
 1400 continue
 1500 continue
 1600 do 1603 i = 1,29
         ti = dfloat(i)/2.9d1
         r2 = x(1)
         t= 1.d0
         do 1601 k = 2,p
              t = t*ti
              r2 = r2 + t*x(k)
 1601    continue
         r2 = -2.d0*r2
         j(i,1) = r2
         t = 1.d0
         r2 = ti*r2
         do 1602 k = 2,p
              j(i,k) = t*(dfloat(k-1) + r2)
              t = t*ti
 1602    continue
 1603 continue
      do 1604 i = 30,31
         do 1604 k = 2,p
 1604         j(i,k) = 0.d0
      j(30,1) = 1.d0
      j(31,1) = -2.d0*x(1)
      j(31,2) = 1.d0
      go to 9999
c  ***  chebyquad  ***
 1700 do 1701 k = 1,n
         tim1 = -1.d0/dfloat(n)
         z = 2.d0*x(k) - 1.d0
         ti = z*tim1
         tpim1 = 0.d0
         tpi = 2.d0*tim1
         z = z + z
         do 1701 i = 1,n
              j(i,k) = tpi
              tpip1 = 4.d0*ti + z*tpi - tpim1
              tpim1 = tpi
              tpi = tpip1
              tip1 = z*ti - tim1
              tim1 = ti
              ti = tip1
 1701         continue
      go to 9999
c  ***  brown and dennis  ***
 1800 do 1801 i = 1, n
         ti = 0.2d0*dfloat(i)
         j(i,1) = 2.0d0*(x(1) + x(2)*ti - dexp(ti))
         j(i,2) = ti*j(i,1)
         t = dsin(ti)
         j(i,3) = 2.0d0*(x(3) + x(4)*t - dcos(ti))
         j(i,4) = t*j(i,3)
 1801    continue
      go to 9999
c  ***  bard  ***
 1900 do 1901 i = 1,15
         j(i,1) = -1.d0
         u = dfloat(i)
         v = 1.6d1 - u
         w = dmin1 (u,v)
         t = u/(x(2)*v + x(3)*w)**2
         j(i,2) = v*t
         j(i,3) = w*t
 1901 continue
      go to 9999
c  *** jennrich & sampson  ***
 2000 do 2001 i = 1,10
         ti = dfloat(i)
         j(i,1) = -ti*dexp(ti*x(1))
         j(i,2) = -ti*dexp(ti*x(2))
 2001    continue
      go to 9999
c  ***  kowalik & osborne  ***
 2100 do 2101 i = 1,11
         t = -1.d0/(ukow(i)**2 + x(3)*ukow(i) + x(4))
         j(i,1) = t*(ukow(i)**2 + x(2)*ukow(i))
         j(i,2) = x(1)*ukow(i)*t
         t = t*j(i,1)*x(1)
         j(i,3) = ukow(i)*t
         j(i,4) = t
 2101 continue
      go to 9999
c  ***  osborne 1  ***
 2200 do 2201 i = 1,33
         ti = 1.0d1*dfloat(1-i)
         j(i,1) = -1.d0
         j(i,2) = -dexp(x(4)*ti)
         j(i,3) = -dexp(x(5)*ti)
         j(i,4) = ti*x(2)*j(i,2)
         j(i,5) = ti*x(3)*j(i,3)
 2201    continue
      go to 9999
c  ***  osborne 2  ***
c     ***  uftolg is a machine-dependent constant.  it is just slightly
c     ***  larger than the log of the smallest positive machine number.
 2300 if (uftolg .eq. 0.d0) uftolg = 1.999d0 * dlog(rmdcon(2))
      do 2302 i = 1,65
         ti = dfloat(1 - i)*1.d-1
         j(i,1) = -dexp(x(5)*ti)
         j(i,5) = x(1)*ti*j(i,1)
         do 2301 k = 2,4
              t = x(k + 7) + ti
              r2 = 0.d0
              theta = -x(k+4)*t*t
              if (theta .gt. uftolg) r2 = -dexp(theta)
              j(i,k) = r2
              r2 = -t*r2*x(k)
              j(i,k+4) = r2*t
              j(i,k+7) = 2.d0*x(k+4)*r2
 2301         continue
 2302    continue
      go to 9999
c  ***  madsen  ***
 2400 j(1,1) = 2.d0*x(1) + x(2)
      j(1,2) = 2.d0*x(2) + x(1)
      j(2,1) = dcos(x(1))
      j(2,2) = 0.d0
      j(3,1) = 0.d0
      j(3,2) = -dsin(x(2))
      go to 9999
c  ***  meyer  ***
 2500 do 2501 i = 1, 16
         ti = dfloat(5*i + 45)
         u = ti + x(3)
         t = dexp(x(2)/u)
         j(i,1) = t
         j(i,2) = x(1)*t/u
         j(i,3) = -x(1)*x(2)*t/(u*u)
 2501    continue
      go to 9999
c  ***  brown  ***
 2600 continue
 2700 continue
 2800 continue
 2900 nm1 = n - 1
      do 2901 k = 1, n
         do 2901 i = 1, nm1
              j(i,k) = 1.0d0
              if (i .eq. k) j(i,k) = 2.0d0
 2901         continue
      do 2903 k = 1, n
         t = 1.0d0
         do 2902 i = 1,n
              if (i .ne. k) t = t*x(i)
 2902         continue
         j(n,k) = t
 2903    continue
      go to 9999
c
c
 9999 return
      end
