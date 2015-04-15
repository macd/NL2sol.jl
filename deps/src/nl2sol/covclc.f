      subroutine covclc(covirc, d, iv, j, n, nn, p, r, v, x)
c
c  ***  compute covariance matrix for nl2itr (nl2sol version 2.2)  ***
c
c  ***  let k = iabs(iv(covreq).  for k .le. 2, a finite-difference
c  ***  hessian h is computed (using func. and grad. values if
c  ***  iv(covreq) is nonnegative, and using only func. values if
c  ***  iv(covreq) is negative).  for scale = 2*f(x) / max(1, n-p),
c  ***  where 2*f(x) is the residual sum of squares, covclc computes...
c  ***             k = 0 or 1...  scale * h**-1 * (j**t * j) * h**-1.
c  ***             k = 2...  scale * h**-1.
c  ***             k .ge. 3...  scale * (j**t * j)**-1.
c
c  ***  parameter declarations  ***
c
      integer covirc, iv(1), n, nn, p
      double precision d(p), j(nn,p), r(n), v(1), x(p)
c     dimension iv(*), v(*)
c
c  ***  local variables  ***
c
      logical havej
      integer cov, gp, gsave1, g1, hc, hmi, hpi, hpm, i, ipivi, ipivk,
     1        ip1, irc, k, kind, kl, l, m, mm1, mm1o2, pp1o2, qtr1,
     2        rd1, stpi, stpm, stp0, wl, w0, w1
      double precision del, half, negpt5, one, t, two, wk, zero
c
c  ***  intrinsic functions  ***
c/+
      integer iabs, max0
      real float
      double precision dabs, dmax1
c/
c  ***  external subroutines  ***
c
      external linvrt, litvmu, livmul, lsqrt, ltsqar, qrfact,
     1         vcopy, vscopy
c
c linvrt... invert lower triangular matrix.
c litvmu... apply inverse-transpose of compact lower triang. matrix.
c livmul... apply inverse of compact lower triang. matrix.
c lsqrt.... compute cholesky factor of (lower trinag. of) a sym. matrix.
c ltsqar... given lower triang. matrix l, compute (l**t)*l.
c qrfact... compute qr decomposition of a matrix.
c vcopy.... copy one vector to another.
c vscopy... set all elements of a vector to a scalar.
c
c  ***  subscripts for iv and v  ***
c
      integer covmat, covreq, delta, delta0, dltfdc, f, fx, g, h, ierr,
     1        ipivot, ipiv0, kagqt, kalm, lmat, mode, nfgcal, qtr,
     2        rd, rsave, savei, switch, toobig, w, xmsave
c
c/6
c      data half/0.5d+0/, negpt5/-0.5d+0/, one/1.d+0/, two/2.d+0/,
c     1     zero/0.d+0/
c/7
      parameter (half=0.5d+0, negpt5=-0.5d+0, one=1.d+0, two=2.d+0,
     1     zero=0.d+0)
c/
c
c/6
c      data covmat/26/, covreq/15/, delta/50/, delta0/44/,
c     1     dltfdc/40/, f/10/, fx/46/, g/28/, h/44/, ierr/32/,
c     2     ipivot/61/, ipiv0/60/, kagqt/35/, kalm/36/,
c     3     lmat/58/, mode/38/, nfgcal/7/, qtr/49/,
c     4     rd/51/, rsave/52/, savei/54/, switch/12/,
c     5     toobig/2/, w/59/, xmsave/49/
c/7
      parameter (covmat=26, covreq=15, delta=50, delta0=44,
     1     dltfdc=40, f=10, fx=46, g=28, h=44, ierr=32,
     2     ipivot=61, ipiv0=60, kagqt=35, kalm=36,
     3     lmat=58, mode=38, nfgcal=7, qtr=49,
     4     rd=51, rsave=52, savei=54, switch=12,
     5     toobig=2, w=59, xmsave=49)
c/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
      covirc = 4
      kind = iv(covreq)
      m = iv(mode)
      if (m .gt. 0) go to 10
         iv(kagqt) = -1
         if (iv(kalm) .gt. 0) iv(kalm) = 0
         if (iabs(kind) .ge. 3) go to 300
         v(fx) = v(f)
         k = iv(rsave)
         call vcopy(n, v(k), r)
 10   if (m .gt. p) go to 200
      if (kind .lt. 0) go to 100
c
c  ***  compute finite-difference hessian using both function and
c  ***  gradient values.
c
      gsave1 = iv(w) + p
      g1 = iv(g)
      if (m .gt. 0) go to 15
c        ***  first call on covclc.  set gsave = g, take first step  ***
         call vcopy(p, v(gsave1), v(g1))
         iv(switch) = iv(nfgcal)
         go to 80
c
 15   del = v(delta)
      x(m) = v(xmsave)
      if (iv(toobig) .eq. 0) go to 30
c
c     ***  handle oversize v(delta)  ***
c
         if (del*x(m) .gt. zero) go to 20
c             ***  we already tried shrinking v(delta), so quit  ***
              iv(covmat) = -2
              go to 190
c
c        ***  try shrinking v(delta)  ***
 20      del = negpt5 * del
         go to 90
c
 30   cov = iv(lmat)
      gp = g1 + p - 1
c
c  ***  set  g = (g - gsave)/del  ***
c
      do 40 i = g1, gp
         v(i) = (v(i) - v(gsave1)) / del
         gsave1 = gsave1 + 1
 40      continue
c
c  ***  add g as new col. to finite-diff. hessian matrix  ***
c
      k = cov + m*(m-1)/2
      l = k + m - 2
      if ( m .eq. 1) go to 60
c
c  ***  set  h(i,m) = 0.5 * (h(i,m) + g(i))  for i = 1 to m-1  ***
c
      do 50 i = k, l
         v(i) = half * (v(i) + v(g1))
         g1 = g1 + 1
 50      continue
c
c  ***  add  h(i,m) = g(i)  for i = m to p  ***
c
 60   l = l + 1
      do 70 i = m, p
         v(l) = v(g1)
         l = l + i
         g1 = g1 + 1
 70      continue
c
 80   m = m + 1
      iv(mode) = m
      if (m .gt. p) go to 190
c
c  ***  choose next finite-difference step, return to get g there  ***
c
      del = v(delta0) * dmax1(one/d(m), dabs(x(m)))
      if (x(m) .lt. zero) del = -del
      v(xmsave) = x(m)
 90   x(m) = x(m) + del
      v(delta) = del
      covirc = 2
      go to 999
c
c  ***  compute finite-difference hessian using function values only.
c
 100  stp0 = iv(w) + p - 1
      mm1 = m - 1
      mm1o2 = m*mm1/2
      if (m .gt. 0) go to 105
c        ***  first call on covclc.  ***
         iv(savei) = 0
         go to 180
c
 105  i = iv(savei)
      if (i .gt. 0) go to 160
      if (iv(toobig) .eq. 0) go to 120
c
c     ***  handle oversize step  ***
c
         stpm = stp0 + m
         del = v(stpm)
         if (del*x(xmsave) .gt. zero) go to 110
c             ***  we already tried shrinking the step, so quit  ***
              iv(covmat) = -2
              go to 999
c
c        ***  try shrinking the step  ***
 110     del = negpt5 * del
         x(m) = x(xmsave) + del
         v(stpm) = del
         covirc = 1
         go to 999
c
c  ***  save f(x + stp(m)*e(m)) in h(p,m)  ***
c
 120  pp1o2 = p * (p-1) / 2
      cov = iv(lmat)
      hpm = cov + pp1o2 + mm1
      v(hpm) = v(f)
c
c  ***  start computing row m of the finite-difference hessian h.  ***
c
      hmi = cov + mm1o2
      if (mm1 .eq. 0) go to 140
      hpi = cov + pp1o2
      do 130 i = 1, mm1
         v(hmi) = v(fx) - (v(f) + v(hpi))
         hmi = hmi + 1
         hpi = hpi + 1
 130     continue
 140  v(hmi) = v(f) - two*v(fx)
c
c  ***  compute function values needed to complete row m of h.  ***
c
      i = 1
c
 150  iv(savei) = i
      stpi = stp0 + i
      v(delta) = x(i)
      x(i) = x(i) + v(stpi)
      if (i .eq. m) x(i) = v(xmsave) - v(stpi)
      covirc = 1
      go to 999
c
 160  x(i) = v(delta)
      if (iv(toobig) .eq. 0) go to 170
c        ***  punt in the event of an oversize step  ***
         iv(covmat) = -2
         go to 999
c
c  ***  finish computing h(m,i)  ***
c
 170  stpi = stp0 + i
      hmi = cov + mm1o2 + i - 1
      stpm = stp0 + m
      v(hmi) = (v(hmi) + v(f)) / (v(stpi)*v(stpm))
      i = i + 1
      if (i .le. m) go to 150
      iv(savei) = 0
      x(m) = v(xmsave)
c
 180  m = m + 1
      iv(mode) = m
      if (m .gt. p) go to 190
c
c  ***  prepare to compute row m of the finite-difference hessian h.
c  ***  compute m-th step size stp(m), then return to obtain
c  ***  f(x + stp(m)*e(m)), where e(m) = m-th std. unit vector.
c
      del = v(dltfdc) * dmax1(one/d(m), dabs(x(m)))
      if (x(m) .lt. zero) del = -del
      v(xmsave) = x(m)
      x(m) = x(m) + del
      stpm = stp0 + m
      v(stpm) = del
      covirc = 1
      go to 999
c
c  ***  restore r, v(f), etc.  ***
c
 190  k = iv(rsave)
      call vcopy(n, r, v(k))
      v(f) = v(fx)
      if (kind .lt. 0) go to 200
         iv(nfgcal) = iv(switch)
         qtr1 = iv(qtr)
         call vcopy(n, v(qtr1), r)
         if (iv(covmat) .lt. 0) go to 999
         covirc = 3
         go to 999
c
 200  cov = iv(lmat)
c
c  ***  the complete finite-diff. hessian is now stored at v(cov).   ***
c  ***  use it to compute the requested covariance matrix.           ***
c
c     ***  compute cholesky factor c of h = c*(c**t)  ***
c     ***  and store it at v(hc).  ***
c
      hc = cov
      if (iabs(kind) .eq. 2) go to 210
         hc = iabs(iv(h))
         iv(h) = -hc
 210  call lsqrt(1, p, v(hc), v(cov), irc)
      iv(covmat) = -1
      if (irc .ne. 0) go to 999
c
      w1 = iv(w) + p
      if (iabs(kind) .gt. 1) go to 350
c
c  ***  covariance = scale * h**-1 * (j**t * j) * h**-1  ***
c
      call vscopy(p*(p+1)/2, v(cov), zero)
      havej = iv(kalm) .eq. (-1)
c     ***  havej = .true. means j is in its original form, while
c     ***  havej = .false. means qrfact has been applied to j.
c
      m = p
      if (havej) m = n
      w0 = w1 - 1
      rd1 = iv(rd)
      do 290 i = 1, m
         if (havej) go to 240
c
c        ***  set w = ipivot * (row i of r matrix from qrfact).  ***
c
              call vscopy(p, v(w1), zero)
              ipivi = ipiv0 + i
              l = w0 + iv(ipivi)
              v(l) = v(rd1)
              rd1 = rd1 + 1
              if (i .eq. p) go to 260
              ip1 = i + 1
              do 230 k = ip1, p
                   ipivk = ipiv0 + k
                   l = w0 + iv(ipivk)
                   v(l) = j(i,k)
 230               continue
              go to 260
c
c        ***  set w = (row i of j).  ***
c
 240     l = w0
         do 250 k = 1, p
              l = l + 1
              v(l) = j(i,k)
 250          continue
c
c        ***  set w = h**-1 * w.  ***
c
 260     call livmul(p, v(w1), v(hc), v(w1))
         call litvmu(p, v(w1), v(hc), v(w1))
c
c        ***  add  w * w**t  to covariance matrix.  ***
c
         kl = cov
         do 280 k = 1, p
              l = w0 + k
              wk = v(l)
              do 270 l = 1, k
                   wl = w0 + l
                   v(kl) = v(kl)  +  wk * v(wl)
                   kl = kl + 1
 270               continue
 280          continue
 290     continue
      go to 380
c
c  ***  covariance = scale * (j**t * j)**-1.  ***
c
 300  rd1 = iv(rd)
      if (iv(kalm) .ne. (-1)) go to 310
c
c        ***  apply qrfact to j  ***
c
         qtr1 = iv(qtr)
         call vcopy(n, v(qtr1), r)
         w1 = iv(w) + p
         call qrfact(nn, n, p, j, v(rd1), iv(ipivot), iv(ierr), 0,
     1               v(w1))
         iv(kalm) = -2
 310  iv(covmat) = -1
      if (iv(ierr) .ne. 0) go to 999
      cov = iv(lmat)
      hc = iabs(iv(h))
      iv(h) = -hc
c
c     ***  set hc = (r matrix from qrfact).  ***
c
      l = hc
      do 340 i = 1, p
         if (i .gt. 1) call vcopy(i-1, v(l), j(1,i))
         l = l + i - 1
         v(l) = v(rd1)
         l = l + 1
         rd1 = rd1 + 1
 340     continue
c
c  ***  the cholesky factor c of the unscaled inverse covariance matrix
c  ***  (or permutation thereof) is stored at v(hc).
c
c  ***  set c = c**-1.
c
 350  call linvrt(p, v(hc), v(hc))
c
c  ***  set c = c**t * c.
c
      call ltsqar(p, v(hc), v(hc))
c
      if (hc .eq. cov) go to 380
c
c     ***  c = permuted, unscaled covariance.
c     ***  set cov = ipivot * c * ipivot**t.
c
         do 370 i = 1, p
              m = ipiv0 + i
              ipivi = iv(m)
              kl = cov-1 + ipivi*(ipivi-1)/2
              do 360 k = 1, i
                   m = ipiv0 + k
                   ipivk = iv(m)
                   l = kl + ipivk
                   if (ipivk .gt. ipivi)
     1                       l = l + (ipivk-ipivi)*(ipivk+ipivi-3)/2
                   v(l) = v(hc)
                   hc = hc + 1
 360               continue
 370          continue
c
 380  iv(covmat) = cov
c
c  ***  apply scale factor = (resid. sum of squares) / max(1,n-p).
c
      t = v(f) / (half * float(max0(1,n-p)))
      k = cov - 1 + p*(p+1)/2
      do 390 i = cov, k
 390     v(i) = t * v(i)
c
 999  return
c  ***  last card of covclc follows  ***
      end
