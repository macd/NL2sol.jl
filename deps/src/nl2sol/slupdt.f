      subroutine slupdt(a, cosmin, p, size, step, u, w, wchmtd, wscale,
     1                  y)
c
c  ***  update symmetric  a  so that  a * step = y  ***
c  ***  (lower triangle of  a  stored rowwise       ***
c
c  ***  parameter declarations  ***
c
      integer p
      double precision a(1), cosmin, size, step(p), u(p), w(p),
     1                 wchmtd(p), wscale, y(p)
c     dimension a(p*(p+1)/2)
c
c  ***  local variables  ***
c
      integer i, j, k
      double precision denmin, sdotwm, t, ui, wi
c
c     ***  constants  ***
      double precision half, one, zero
c
c  ***  intrinsic functions  ***
c/+
      double precision dabs, dmin1
c/
c  ***  external functions and subroutines  ***
c
      external dotprd, slvmul, v2norm
      double precision dotprd, v2norm
c
c/6
c      data half/0.5d+0/, one/1.d+0/, zero/0.d+0/
c/7
      parameter (half=0.5d+0, one=1.d+0, zero=0.d+0)
c/
c
c-----------------------------------------------------------------------
c
      sdotwm = dotprd(p, step, wchmtd)
      denmin = cosmin * v2norm(p,step) * v2norm(p,wchmtd)
      wscale = one
      if (denmin .ne. zero) wscale = dmin1(one, dabs(sdotwm/denmin))
      t = zero
      if (sdotwm .ne. zero) t = wscale / sdotwm
      do 10 i = 1, p
 10      w(i) = t * wchmtd(i)
      call slvmul(p, u, a, step)
      t = half * (size * dotprd(p, step, u)  -  dotprd(p, step, y))
      do 20 i = 1, p
 20      u(i) = t*w(i) + y(i) - size*u(i)
c
c  ***  set  a = a + u*(w**t) + w*(u**t)  ***
c
      k = 1
      do 40 i = 1, p
         ui = u(i)
         wi = w(i)
         do 30 j = 1, i
              a(k) = size*a(k) + ui*w(j) + wi*u(j)
              k = k + 1
 30           continue
 40      continue
c
 999  return
c  ***  last card of slupdt follows  ***
      end
