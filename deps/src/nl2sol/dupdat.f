      subroutine dupdat(d, iv, j, n, nn, p, v)
c
c  ***  update scale vector d for nl2itr (nl2sol version 2.2)  ***
c
c  ***  parameter declarations  ***
c
      integer iv(1), n, nn, p
      double precision d(p), j(nn,p), v(1)
c     dimension iv(*), v(*)
c
c  ***  local variables  ***
c
      integer d0, i, jtoli, s1
      double precision sii, t, vdfac
c
c     ***  constants  ***
      double precision zero
c
c  ***  intrinsic functions  ***
c/+
      double precision dmax1, dsqrt
c/
c  ***  external function  ***
c
      external v2norm
      double precision v2norm
c
c  ***  subscripts for iv and v  ***
c
      integer dfac, dtype, jtol0, niter, s
c/6
c      data dfac/41/, dtype/16/, jtol0/86/, niter/31/, s/53/
c/7
      parameter (dfac=41, dtype=16, jtol0=86, niter=31, s=53)
c/
c
c/6
c      data zero/0.d+0/
c/7
      parameter (zero=0.d+0)
c/
c
c-----------------------------------------------------------------------
c
      i = iv(dtype)
      if (i .eq. 1) go to 20
         if (iv(niter) .gt. 0) go to 999
c
 20   vdfac = v(dfac)
      d0 = jtol0 + p
      s1 = iv(s) - 1
      do 30 i = 1, p
         s1 = s1 + i
         sii = v(s1)
         t = v2norm(n, j(1,i))
         if (sii .gt. zero) t = dsqrt(t*t + sii)
         jtoli = jtol0 + i
         d0 = d0 + 1
         if (t .lt. v(jtoli)) t = dmax1(v(d0), v(jtoli))
         d(i) = dmax1(vdfac*d(i), t)
 30      continue
c
 999  return
c  ***  last card of dupdat follows  ***
      end
