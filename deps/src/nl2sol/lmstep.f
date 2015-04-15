      subroutine lmstep(d, g, ierr, ipivot, ka, p, qtr, r, step, v, w)
c
c  ***  compute levenberg-marquardt step using more-hebden technique  **
c  ***  nl2sol version 2.2.  ***
c
c  ***  parameter declarations  ***
c
      integer ierr, ka, p
      integer ipivot(p)
      double precision d(p), g(p), qtr(p), r(1), step(p), v(21), w(1)
c     dimension w(p*(p+5)/2 + 4)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  purpose  ***
c
c        given the r matrix from the qr decomposition of a jacobian
c     matrix, j, as well as q-transpose times the corresponding
c     residual vector, resid, this subroutine computes a levenberg-
c     marquardt step of approximate length v(radius) by the more-
c     technique.
c
c  ***  parameter description  ***
c
c      d (in)  = the scale vector.
c      g (in)  = the gradient vector (j**t)*r.
c   ierr (i/o) = return code from qrfact or qrfgs -- 0 means r has
c             full rank.
c ipivot (i/o) = permutation array from qrfact or qrfgs, which compute
c             qr decompositions with column pivoting.
c     ka (i/o).  ka .lt. 0 on input means this is the first call on
c             lmstep for the current r and qtr.  on output ka con-
c             tains the number of hebden iterations needed to determine
c             step.  ka = 0 means a gauss-newton step.
c      p (in)  = number of parameters.
c    qtr (in)  = (q**t)*resid = q-transpose times the residual vector.
c      r (in)  = the r matrix, stored compactly by columns.
c   step (out) = the levenberg-marquardt step computed.
c      v (i/o) contains various constants and variables described below.
c      w (i/o) = workspace of length p*(p+5)/2 + 4.
c
c  ***  entries in v  ***
c
c v(dgnorm) (i/o) = 2-norm of (d**-1)*g.
c v(dstnrm) (i/o) = 2-norm of d*step.
c v(dst0)   (i/o) = 2-norm of gauss-newton step (for nonsing. j).
c v(epslon) (in) = max. rel. error allowed in twonorm(r)**2 minus
c             twonorm(r - j*step)**2.  (see algorithm notes below.)
c v(gtstep) (out) = inner product between g and step.
c v(nreduc) (out) = half the reduction in the sum of squares predicted
c             for a gauss-newton step.
c v(phmnfc) (in)  = tol. (together with v(phmxfc)) for accepting step
c             (more*s sigma).  the error v(dstnrm) - v(radius) must lie
c             between v(phmnfc)*v(radius) and v(phmxfc)*v(radius).
c v(phmxfc) (in)  (see v(phmnfc).)
c v(preduc) (out) = half the reduction in the sum of squares predicted
c             by the step returned.
c v(radius) (in)  = radius of current (scaled) trust region.
c v(rad0)   (i/o) = value of v(radius) from previous call.
c v(stppar) (i/o) = marquardt parameter (or its negative if the special
c             case mentioned below in the algorithm notes occurs).
c
c note -- see data statement below for values of above subscripts.
c
c  ***  usage notes  ***
c
c     if it is desired to recompute step using a different value of
c     v(radius), then this routine may be restarted by calling it
c     with all parameters unchanged except v(radius).  (this explains
c     why many parameters are listed as i/o).  on an intiial call (one
c     with ka = -1), the caller need only have initialized d, g, ka, p,
c     qtr, r, v(epslon), v(phmnfc), v(phmxfc), v(radius), and v(rad0).
c
c  ***  application and usage restrictions  ***
c
c     this routine is called as part of the nl2sol (nonlinear least-
c     squares) package (ref. 1).
c
c  ***  algorithm notes  ***
c
c     this code implements the step computation scheme described in
c     refs. 2 and 4.  fast givens transformations (see ref. 3, pp. 60-
c     62) are used to compute step with a nonzero marquardt parameter.
c        a special case occurs if j is (nearly) singular and v(radius)
c     is sufficiently large.  in this case the step returned is such
c     that  twonorm(r)**2 - twonorm(r - j*step)**2  differs from its
c     optimal value by less than v(epslon) times this optimal value,
c     where j and r denote the original jacobian and residual.  (see
c     ref. 2 for more details.)
c
c  ***  functions and subroutines called  ***
c
c dotprd - returns inner product of two vectors.
c litvmu - apply inverse-transpose of compact lower triang. matrix.
c livmul - apply inverse of compact lower triang. matrix.
c vcopy  - copies one vector to another.
c v2norm - returns 2-norm of a vector.
c
c  ***  references  ***
c
c 1.  dennis, j.e., gay, d.m., and welsch, r.e. (1981), an adaptive
c             nonlinear least-squares algorithm, acm trans. math.
c             software, vol. 7, no. 3.
c 2.  gay, d.m. (1981), computing optimal locally constrained steps,
c             siam j. sci. statist. computing, vol. 2, no. 2, pp.
c             186-197.
c 3.  lawson, c.l., and hanson, r.j. (1974), solving least squares
c             problems, prentice-hall, englewood cliffs, n.j.
c 4.  more, j.j. (1978), the levenberg-marquardt algorithm, implemen-
c             tation and theory, pp.105-116 of springer lecture notes
c             in mathematics no. 630, edited by g.a. watson, springer-
c             verlag, berlin and new york.
c
c  ***  general  ***
c
c     coded by david m. gay.
c     this subroutine was written in connection with research
c     supported by the national science foundation under grants
c     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and
c     mcs-7906671.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  local variables  ***
c
      integer dstsav, i, ip1, i1, j1, k, kalim, l, lk0, phipin,
     1        pp1o2, res, res0, rmat, rmat0, uk0
      double precision a, adi, alphak, b, dfacsq, dst, dtol, d1, d2,
     1                 lk, oldphi, phi, phimax, phimin, psifac, rad,
     2                 si, sj, sqrtak, t, twopsi, uk, wl
c
c     ***  constants  ***
      double precision dfac, eight, half, negone, one, p001, three,
     1                 ttol, zero
c
c  ***  intrinsic functions  ***
c/+
      integer iabs
      double precision dabs, dmax1, dmin1, dsqrt
c/
c  ***  external functions and subroutines  ***
c
      external dotprd, litvmu, livmul, vcopy, v2norm
      double precision dotprd, v2norm
c
c  ***  subscripts for v  ***
c
      integer dgnorm, dstnrm, dst0, epslon, gtstep, nreduc, phmnfc,
     1        phmxfc, preduc, radius, rad0, stppar
c/6
c      data dgnorm/1/, dstnrm/2/, dst0/3/, epslon/19/,
c     1     gtstep/4/, nreduc/6/, phmnfc/20/,
c     2     phmxfc/21/, preduc/7/, radius/8/,
c     3     rad0/9/, stppar/5/
c/7
      parameter (dgnorm=1, dstnrm=2, dst0=3, epslon=19,
     1     gtstep=4, nreduc=6, phmnfc=20,
     2     phmxfc=21, preduc=7, radius=8,
     3     rad0=9, stppar=5)
c/
c
c/6
c      data dfac/256.d+0/, eight/8.d+0/, half/0.5d+0/, negone/-1.d+0/,
c     1     one/1.d+0/, p001/1.d-3/, three/3.d+0/, ttol/2.5d+0/,
c     2     zero/0.d+0/
c/7
      parameter (dfac=256.d+0, eight=8.d+0, half=0.5d+0, negone=-1.d+0,
     1     one=1.d+0, p001=1.d-3, three=3.d+0, ttol=2.5d+0,
     2     zero=0.d+0)
c/
c
c  ***  body  ***
c
c     ***  for use in recomputing step, the final values of lk and uk,
c     ***  the inverse derivative of more*s phi at 0 (for nonsing. j)
c     ***  and the value returned as v(dstnrm) are stored at w(lk0),
c     ***  w(uk0), w(phipin), and w(dstsav) respectively.
      lk0 = p + 1
      phipin = lk0 + 1
      uk0 = phipin + 1
      dstsav = uk0 + 1
      rmat0 = dstsav
c     ***  a copy of the r-matrix from the qr decomposition of j is
c     ***  stored in w starting at w(rmat), and a copy of the residual
c     ***  vector is stored in w starting at w(res).  the loops below
c     ***  that update the qr decomp. for a nonzero marquardt parameter
c     ***  work on these copies.
      rmat = rmat0 + 1
      pp1o2 = p * (p + 1) / 2
      res0 = pp1o2 + rmat0
      res = res0 + 1
      rad = v(radius)
      if (rad .gt. zero)
     1   psifac = v(epslon)/((eight*(v(phmnfc) + one) + three) * rad**2)
      phimax = v(phmxfc) * rad
      phimin = v(phmnfc) * rad
c     ***  dtol, dfac, and dfacsq are used in rescaling the fast givens
c     ***  representation of the updated qr decomposition.
      dtol = one/dfac
      dfacsq = dfac*dfac
c     ***  oldphi is used to detect limits of numerical accuracy.  if
c     ***  we recompute step and it does not change, then we accept it.
      oldphi = zero
      lk = zero
      uk = zero
      kalim = ka + 12
c
c  ***  start or restart, depending on ka  ***
c
      if (ka) 10, 20, 370
c
c  ***  fresh start -- compute v(nreduc)  ***
c
 10   ka = 0
      kalim = 12
      k = p
      if (ierr .ne. 0) k = iabs(ierr) - 1
      v(nreduc) = half*dotprd(k, qtr, qtr)
c
c  ***  set up to try initial gauss-newton step  ***
c
 20   v(dst0) = negone
      if (ierr .ne. 0) go to 90
c
c  ***  compute gauss-newton step  ***
c
c     ***  note -- the r-matrix is stored compactly by columns in
c     ***  r(1), r(2), r(3), ...  it is the transpose of a
c     ***  lower triangular matrix stored compactly by rows, and we
c     ***  treat it as such when using litvmu and livmul.
      call litvmu(p, w, r, qtr)
c     ***  temporarily store permuted -d*step in step.
      do 60 i = 1, p
         j1 = ipivot(i)
         step(i) = d(j1)*w(i)
 60      continue
      dst = v2norm(p, step)
      v(dst0) = dst
      phi = dst - rad
      if (phi .le. phimax) go to 410
c     ***  if this is a restart, go to 110  ***
      if (ka .gt. 0) go to 110
c
c  ***  gauss-newton step was unacceptable.  compute l0  ***
c
      do 70 i = 1, p
         j1 = ipivot(i)
         step(i) = d(j1)*(step(i)/dst)
 70      continue
      call livmul(p, step, r, step)
      t = one / v2norm(p, step)
      w(phipin) = (t/dst)*t
      lk = phi*w(phipin)
c
c  ***  compute u0  ***
c
 90   do 100 i = 1, p
 100     w(i) = g(i)/d(i)
      v(dgnorm) = v2norm(p, w)
      uk = v(dgnorm)/rad
      if (uk .le. zero) go to 390
c
c     ***  alphak will be used as the current marquardt parameter.  we
c     ***  use more*s scheme for initializing it.
      alphak = dabs(v(stppar)) * v(rad0)/rad
c
c
c  ***  top of loop -- increment ka, copy r to rmat, qtr to res  ***
c
 110  ka = ka + 1
      call vcopy(pp1o2, w(rmat), r)
      call vcopy(p, w(res), qtr)
c
c  ***  safeguard alphak and initialize fast givens scale vector.  ***
c
      if (alphak .le. zero .or. alphak .lt. lk .or. alphak .ge. uk)
     1             alphak = uk * dmax1(p001, dsqrt(lk/uk))
      sqrtak = dsqrt(alphak)
      do 120 i = 1, p
 120     w(i) = one
c
c  ***  add alphak*d and update qr decomp. using fast givens trans.  ***
c
      do 270 i = 1, p
c        ***  generate, apply 1st givens trans. for row i of alphak*d.
c        ***  (use step to store temporary row)  ***
         l = i*(i+1)/2 + rmat0
         wl = w(l)
         d2 = one
         d1 = w(i)
         j1 = ipivot(i)
         adi = sqrtak*d(j1)
         if (adi .ge. dabs(wl)) go to 150
 130     a = adi/wl
         b = d2*a/d1
         t = a*b + one
         if (t .gt. ttol) go to 150
         w(i) = d1/t
         d2 = d2/t
         w(l) = t*wl
         a = -a
         do 140 j1 = i, p
              l = l + j1
              step(j1) = a*w(l)
 140          continue
         go to 170
c
 150     b = wl/adi
         a = d1*b/d2
         t = a*b + one
         if (t .gt. ttol) go to 130
         w(i) = d2/t
         d2 = d1/t
         w(l) = t*adi
         do 160 j1 = i, p
              l = l + j1
              wl = w(l)
              step(j1) = -wl
              w(l) = a*wl
 160          continue
c
 170     if (i .eq. p) go to 280
c
c        ***  now use givens trans. to zero elements of temp. row  ***
c
         ip1 = i + 1
         do 260 i1 = ip1, p
              l = i1*(i1+1)/2 + rmat0
              wl = w(l)
              si = step(i1-1)
              d1 = w(i1)
c
c             ***  rescale row i1 if necessary  ***
c
              if (d1 .ge. dtol) go to 190
                   d1 = d1*dfacsq
                   wl = wl/dfac
                   k = l
                   do 180 j1 = i1, p
                        k = k + j1
                        w(k) = w(k)/dfac
 180                    continue
c
c             ***  use givens trans. to zero next element of temp. row
c
 190          if (dabs(si) .gt. dabs(wl)) go to 220
              if (si .eq. zero) go to 260
 200          a = si/wl
              b = d2*a/d1
              t = a*b + one
              if (t .gt. ttol) go to 220
              w(l) = t*wl
              w(i1) = d1/t
              d2 = d2/t
              do 210 j1 = i1, p
                   l = l + j1
                   wl = w(l)
                   sj = step(j1)
                   w(l) = wl + b*sj
                   step(j1) = sj - a*wl
 210               continue
              go to 240
c
 220          b = wl/si
              a = d1*b/d2
              t = a*b + one
              if (t .gt. ttol) go to 200
              w(i1) = d2/t
              d2 = d1/t
              w(l) = t*si
              do 230 j1 = i1, p
                   l = l + j1
                   wl = w(l)
                   sj = step(j1)
                   w(l) = a*wl + sj
                   step(j1) = b*sj - wl
 230               continue
c
c             ***  rescale temp. row if necessary  ***
c
 240          if (d2 .ge. dtol) go to 260
                   d2 = d2*dfacsq
                   do 250 k = i1, p
 250                    step(k) = step(k)/dfac
 260          continue
 270     continue
c
c  ***  compute step  ***
c
 280  call litvmu(p, w(res), w(rmat), w(res))
c     ***  recover step and store permuted -d*step at w(res)  ***
      do 290 i = 1, p
         j1 = ipivot(i)
         k = res0 + i
         t = w(k)
         step(j1) = -t
         w(k) = t*d(j1)
 290     continue
      dst = v2norm(p, w(res))
      phi = dst - rad
      if (phi .le. phimax .and. phi .ge. phimin) go to 430
      if (oldphi .eq. phi) go to 430
      oldphi = phi
c
c  ***  check for (and handle) special case  ***
c
      if (phi .gt. zero) go to 310
         if (ka .ge. kalim) go to 430
              twopsi = alphak*dst*dst - dotprd(p, step, g)
              if (alphak .ge. twopsi*psifac) go to 310
                   v(stppar) = -alphak
                   go to 440
c
c  ***  unacceptable step -- update lk, uk, alphak, and try again  ***
c
 300  if (phi .lt. zero) uk = dmin1(uk, alphak)
      go to 320
 310  if (phi .lt. zero) uk = alphak
 320  do 330 i = 1, p
         j1 = ipivot(i)
         k = res0 + i
         step(i) = d(j1) * (w(k)/dst)
 330     continue
      call livmul(p, step, w(rmat), step)
      do 340 i = 1, p
 340     step(i) = step(i) / dsqrt(w(i))
      t = one / v2norm(p, step)
      alphak = alphak + t*phi*t/rad
      lk = dmax1(lk, alphak)
      go to 110
c
c  ***  restart  ***
c
 370  lk = w(lk0)
      uk = w(uk0)
      if (v(dst0) .gt. zero .and. v(dst0) - rad .le. phimax) go to 20
      alphak = dabs(v(stppar))
      dst = w(dstsav)
      phi = dst - rad
      t = v(dgnorm)/rad
      if (rad .gt. v(rad0)) go to 380
c
c        ***  smaller radius  ***
         uk = t
         if (alphak .le. zero) lk = zero
         if (v(dst0) .gt. zero) lk = dmax1(lk, (v(dst0)-rad)*w(phipin))
         go to 300
c
c     ***  bigger radius  ***
 380  if (alphak .le. zero .or. uk .gt. t) uk = t
      lk = zero
      if (v(dst0) .gt. zero) lk = dmax1(lk, (v(dst0)-rad)*w(phipin))
      go to 300
c
c  ***  special case -- rad .le. 0 or (g = 0 and j is singular)  ***
c
 390  v(stppar) = zero
      dst = zero
      lk = zero
      uk = zero
      v(gtstep) = zero
      v(preduc) = zero
      do 400 i = 1, p
 400     step(i) = zero
      go to 450
c
c  ***  acceptable gauss-newton step -- recover step from w  ***
c
 410  alphak = zero
      do 420 i = 1, p
         j1 = ipivot(i)
         step(j1) = -w(i)
 420     continue
c
c  ***  save values for use in a possible restart  ***
c
 430  v(stppar) = alphak
 440  v(gtstep) = dotprd(p, step, g)
      v(preduc) = half * (alphak*dst*dst - v(gtstep))
 450  v(dstnrm) = dst
      w(dstsav) = dst
      w(lk0) = lk
      w(uk0) = uk
      v(rad0) = rad
c
 999  return
c
c  ***  last card of lmstep follows  ***
      end
