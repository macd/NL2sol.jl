      subroutine nl2itr (d, iv, j, n, nn, p, r, v, x)
c
c  ***  carry out nl2sol (nonlinear least-squares) iterations  ***
c  ***  (nl2sol version 2.2)  ***
c
c  ***  parameter declarations  ***
c
      integer iv(1), n, nn, p
      double precision d(p), j(nn,p), r(n), v(1), x(p)
c     dimension iv(60+p), v(93 + 2*n + p*(3*p+31)/2)
c
c
c--------------------------  parameter usage  --------------------------
c
c d.... scale vector.
c iv... integer value array.
c j.... n by p jacobian matrix (lead dimension nn).
c n.... number of observations (components in r).
c nn... lead dimension of j.
c p.... number of parameters (components in x).
c r.... residual vector.
c v.... floating-point value array.
c x.... parameter vector.
c
c  ***  discussion  ***
c
c        parameters iv, n, p, v, and x are the same as the correspond-
c     ing ones to nl2sol (which see), except that v can be shorter
c     (since the part of v that nl2sol uses for storing d, j, and r is
c     not needed).  moreover, compared with nl2sol, iv(1) may have the
c     two additional output values 1 and 2, which are explained below,
c     as is the use of iv(toobig) and iv(nfgcal).  the values iv(d),
c     iv(j), and iv(r), which are output values from nl2sol (and
c     nl2sno), are not referenced by nl2itr or the subroutines it calls.
c        on a fresh start, i.e., a call on nl2itr with iv(1) = 0 or 12,
c     nl2itr assumes that r = r(x), the residual at x, and j = j(x),
c     the corresponding jacobian matrix of r at x.
c
c iv(1) = 1 means the caller should set r to r(x), the residual at x,
c             and call nl2itr again, having changed none of the other
c             parameters.  an exception occurs if r cannot be evaluated
c             at x (e.g. if r would overflow), which may happen because
c             of an oversized step.  in this case the caller should set
c             iv(toobig) = iv(2) to 1, which will cause nl2itr to ig-
c             nore r and try a smaller step.  the parameter nf that
c             nl2sol passes to calcr (for possible use by calcj) is a
c             copy of iv(nfcall) = iv(6).
c iv(1) = 2 means the caller should set j to j(x), the jacobian matrix
c             of r at x, and call nl2itr again.  the caller may change
c             d at this time, but should not change any of the other
c             parameters.  the parameter nf that nl2sol passes to
c             calcj is iv(nfgcal) = iv(7).  if j cannot be evaluated
c             at x, then the caller may set iv(nfgcal) to 0, in which
c             case nl2itr will return with iv(1) = 15.
c
c  ***  general  ***
c
c     coded by david m. gay.
c     this subroutine was written in connection with research
c     supported by the national science foundation under grants
c
c     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and
c     mcs-7906671.
c        (see nl2sol for references.)
c
c+++++++++++++++++++++++++++  declarations  ++++++++++++++++++++++++++++
c
c  ***  local variables  ***
c
      integer dummy, dig1, g1, g01, h0, h1, i, im1, ipivi, ipivk, ipiv1,
     1        ipk, k, km1, l, lky1, lmat1, lstgst, m, pp1o2, qtr1,
     2        rdk, rd0, rd1, rsave1, smh, sstep, step1, stpmod, s1,
     3        temp1, temp2, w1, x01
      double precision e, rdof1, sttsst, t, t1
c
c     ***  constants  ***
c
      double precision half, negone, one, zero
c
c  ***  intrinsic functions  ***
c/+
      integer iabs
      double precision dabs
c/
c  ***  external functions and subroutines  ***
c
      external assess, covclc, dotprd, dupdat, gqtstp, itsmry, lmstep,
     1         parchk, qapply, qrfact, rptmul, slupdt, slvmul, stopx,
     2         vaxpy, vcopy, vscopy, v2norm
      logical stopx
      double precision dotprd, v2norm
c
c assess... assesses candidate step.
c covclc... computes covariance matrix.
c dotprd... returns inner product of two vectors.
c dupdat... updates scale vector d.
c gqtstp... computes goldfeld-quandt-trotter step (augmented model).
c itsmry... prints iteration summary and info about initial and final x.
c lmstep... computes levenberg-marquardt step (gauss-newton model).
c parchk... checks validity of input iv and v values.
c qapply... applies orthogonal matrix q from qrfact to a vector.
c qrfact... computes qr decomposition of a matrix via householder trans.
c rptmul... multiplies vector by the r matrix (and/or its transpose)
c             stored by qrfact.
c slupdt... performs quasi-newton update on compactly stored lower tri-
c             angle of a symmetric matrix.
c stopx.... returns .true. if the break key has been pressed.
c vaxpy.... computes scalar times one vector plus another.
c vcopy.... copies one vector to another.
c vscopy... sets all elements of a vector to a scalar.
c v2norm... returns the 2-norm of a vector.
c
c  ***  subscripts for iv and v  ***
c
      integer cnvcod, cosmin, covmat, covprt, covreq, dgnorm, dig,
     1        dinit, dstnrm, dtype, d0init, f, fdif, fuzz,
     2        f0, g, gtstep, h, ierr, incfac, inits, ipivot, ipiv0, irc,
     3        jtinit, jtol1, kagqt, kalm, lky, lmat, lmax0, mode, model,
     4        mxfcal, mxiter, nfcall, nfgcal, nfcov, ngcov, ngcall,
     5        niter, nvsave, phmxfc, preduc, qtr, radfac, radinc,
     6        radius, rad0, rd, restor, rlimit, rsave, s, size, step,
     7        stglim, stlstg, stppar, sused, switch, toobig, tuner4,
     8        tuner5, vsave1, w, wscale, xirc, x0
c
c  ***  iv subscript values  ***
c
c/6
c      data cnvcod/34/, covmat/26/, covprt/14/,
c     1     covreq/15/, dig/43/, dtype/16/, g/28/, h/44/,
c     2     ierr/32/, inits/25/, ipivot/61/, ipiv0/60/,
c     3     irc/3/, kagqt/35/, kalm/36/, lky/37/, lmat/58/,
c     4     mode/38/, model/5/, mxfcal/17/, mxiter/18/,
c     5     nfcall/6/, nfgcal/7/, nfcov/40/, ngcov/41/,
c     6     ngcall/30/, niter/31/, qtr/49/,
c     7     radinc/8/, rd/51/, restor/9/, rsave/52/, s/53/,
c     8     step/55/, stglim/11/, stlstg/56/, sused/57/,
c     9     switch/12/, toobig/2/, w/59/, xirc/13/, x0/60/
c/7
      parameter (cnvcod=34, covmat=26, covprt=14,
     1     covreq=15, dig=43, dtype=16, g=28, h=44,
     2     ierr=32, inits=25, ipivot=61, ipiv0=60,
     3     irc=3, kagqt=35, kalm=36, lky=37, lmat=58,
     4     mode=38, model=5, mxfcal=17, mxiter=18,
     5     nfcall=6, nfgcal=7, nfcov=40, ngcov=41,
     6     ngcall=30, niter=31, qtr=49,
     7     radinc=8, rd=51, restor=9, rsave=52, s=53,
     8     step=55, stglim=11, stlstg=56, sused=57,
     9     switch=12, toobig=2, w=59, xirc=13, x0=60)
c/
c
c  ***  v subscript values  ***
c
c/6
c      data cosmin/43/, dgnorm/1/, dinit/38/, dstnrm/2/,
c     1     d0init/37/, f/10/, fdif/11/, fuzz/45/,
c     2     f0/13/, gtstep/4/, incfac/23/,
c     3     jtinit/39/, jtol1/87/, lmax0/35/,
c     4     nvsave/9/, phmxfc/21/, preduc/7/,
c     5     radfac/16/, radius/8/, rad0/9/, rlimit/42/,
c     6     size/47/, stppar/5/, tuner4/29/, tuner5/30/,
c     7     vsave1/78/, wscale/48/
c/7
      parameter (cosmin=43, dgnorm=1, dinit=38, dstnrm=2,
     1     d0init=37, f=10, fdif=11, fuzz=45,
     2     f0=13, gtstep=4, incfac=23,
     3     jtinit=39, jtol1=87, lmax0=35,
     4     nvsave=9, phmxfc=21, preduc=7,
     5     radfac=16, radius=8, rad0=9, rlimit=42,
     6     size=47, stppar=5, tuner4=29, tuner5=30,
     7     vsave1=78, wscale=48)
c/
c
c
c/6
c      data half/0.5d+0/, negone/-1.d+0/, one/1.d+0/, zero/0.d+0/
c/7
      parameter (half=0.5d+0, negone=-1.d+0, one=1.d+0, zero=0.d+0)
c/
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      i = iv(1)
      if (i .eq. 1) go to 20
      if (i .eq. 2) go to 50
c
c  ***  check validity of iv and v input values  ***
c
c     ***  note -- if iv(1) = 0, then parchk calls dfault(iv, v)  ***
      call parchk(iv, n, nn, p, v)
      i = iv(1) - 2
      if (i .gt. 10) go to 999
      go to (350, 350, 350, 350, 350, 350, 195, 160, 195, 10), i
c
c  ***  initialization and storage allocation  ***
c
 10   iv(niter) = 0
      iv(nfcall) = 1
      iv(ngcall) = 1
      iv(nfgcal) = 1
      iv(mode) = -1
      iv(stglim) = 2
      iv(toobig) = 0
      iv(cnvcod) = 0
      iv(covmat) = 0
      iv(nfcov) = 0
      iv(ngcov) = 0
      iv(kalm) = -1
      iv(radinc) = 0
      iv(s) = jtol1 + 2*p
      pp1o2 = p * (p + 1) / 2
      iv(x0) = iv(s) + pp1o2
      iv(step) = iv(x0) + p
      iv(stlstg) = iv(step) + p
      iv(dig) = iv(stlstg) + p
      iv(g) = iv(dig) + p
      iv(lky) = iv(g) + p
      iv(rd) = iv(lky) + p
      iv(rsave) = iv(rd) + p
      iv(qtr) = iv(rsave) + n
      iv(h) = iv(qtr) + n
      iv(w) = iv(h) + pp1o2
      iv(lmat) = iv(w) + 4*p + 7
c     +++ length of w = p*(p+9)/2 + 7.  lmat is contained in w.
      if (v(dinit) .ge. zero) call vscopy(p, d, v(dinit))
      if (v(jtinit) .gt. zero) call vscopy(p, v(jtol1), v(jtinit))
      i = jtol1 + p
      if (v(d0init) .gt. zero) call vscopy(p, v(i), v(d0init))
      v(rad0) = zero
      v(stppar) = zero
      v(radius) = v(lmax0) / (one + v(phmxfc))
c
c  ***  set initial model and s matrix  ***
c
      iv(model) = 1
      if (iv(inits) .eq. 2) iv(model) = 2
      s1 = iv(s)
      if (iv(inits) .eq. 0) call vscopy(pp1o2, v(s1), zero)
c
c  ***  compute function value (half the sum of squares)  ***
c
 20   t = v2norm(n, r)
      if (t .gt. v(rlimit)) iv(toobig) = 1
      if (iv(toobig) .ne. 0) go to 30
      v(f) = half * t**2
 30   if (iv(mode)) 40, 350, 730
c
 40   if (iv(toobig) .eq. 0) go to 60
         iv(1) = 13
         go to 900
c
c  ***  make sure jacobian could be computed  ***
c
 50   if (iv(nfgcal) .ne. 0) go to 60
         iv(1) = 15
         go to 900
c
c  ***  compute gradient  ***
c
 60   iv(kalm) = -1
      g1 = iv(g)
      do 70 i = 1, p
         v(g1) = dotprd(n, r, j(1,i))
         g1 = g1 + 1
 70      continue
      if (iv(mode) .gt. 0) go to 710
c
c  ***  update d and make copies of r for possible use later  ***
c
      if (iv(dtype) .gt. 0) call dupdat(d, iv, j, n, nn, p, v)
      rsave1 = iv(rsave)
      call vcopy(n, v(rsave1), r)
      qtr1 = iv(qtr)
      call vcopy(n, v(qtr1), r)
c
c  ***  compute  d**-1 * gradient  ***
c
      g1 = iv(g)
      dig1 = iv(dig)
      k = dig1
      do 80 i = 1, p
         v(k) = v(g1) / d(i)
         k = k + 1
         g1 = g1 + 1
 80      continue
      v(dgnorm) = v2norm(p, v(dig1))
c
      if (iv(cnvcod) .ne. 0) go to 700
      if (iv(mode) .eq. 0) go to 570
      iv(mode) = 0
c
c
c-----------------------------  main loop  -----------------------------
c
c
c  ***  print iteration summary, check iteration limit  ***
c
 150  call itsmry(d, iv, p, v, x)
 160  k = iv(niter)
      if (k .lt. iv(mxiter)) go to 170
         iv(1) = 10
         go to 900
 170  iv(niter) = k + 1
c
c  ***  update radius  ***
c
      if (k .eq. 0) go to 185
      step1 = iv(step)
      do 180 i = 1, p
         v(step1) = d(i) * v(step1)
         step1 = step1 + 1
 180     continue
      step1 = iv(step)
      v(radius) = v(radfac) * v2norm(p, v(step1))
c
c  ***  initialize for start of next iteration  ***
c
 185  x01 = iv(x0)
      v(f0) = v(f)
      iv(kagqt) = -1
      iv(irc) = 4
      iv(h) = -iabs(iv(h))
      iv(sused) = iv(model)
c
c     ***  copy x to x0  ***
c
      call vcopy(p, v(x01), x)
c
c  ***  check stopx and function evaluation limit  ***
c
 190  if (.not. stopx(dummy)) go to 200
         iv(1) = 11
         go to 205
c
c     ***  come here when restarting after func. eval. limit or stopx.
c
 195  if (v(f) .ge. v(f0)) go to 200
         v(radfac) = one
         k = iv(niter)
         go to 170
c
 200  if (iv(nfcall) .lt. iv(mxfcal) + iv(nfcov)) go to 210
         iv(1) = 9
 205     if (v(f) .ge. v(f0)) go to 900
c
c        ***  in case of stopx or function evaluation limit with
c        ***  improved v(f), evaluate the gradient at x.
c
              iv(cnvcod) = iv(1)
              go to 560
c
c. . . . . . . . . . . . .  compute candidate step  . . . . . . . . . .
c
 210  step1 = iv(step)
      w1 = iv(w)
      if (iv(model) .eq. 2) go to 240
c
c  ***  compute levenberg-marquardt step  ***
c
         qtr1 = iv(qtr)
         if (iv(kalm) .ge. 0) go to 215
              rd1 = iv(rd)
              if (-1 .eq. iv(kalm)) call qrfact(nn, n, p, j, v(rd1),
     1                                   iv(ipivot), iv(ierr), 0, v(w1))
              call qapply(nn, n, p, j, v(qtr1), iv(ierr))
 215     h1 = iv(h)
         if (h1 .gt. 0) go to 230
c
c        ***  copy r matrix to h  ***
c
              h1 = -h1
              iv(h) = h1
              k = h1
              rd1 = iv(rd)
              v(k) = v(rd1)
              if (p .eq. 1) go to 230
              do 220 i = 2, p
                   call vcopy(i-1, v(k+1), j(1,i))
                   k = k + i
                   rd1 = rd1 + 1
                   v(k) = v(rd1)
 220               continue
c
 230     g1 = iv(g)
         call lmstep(d, v(g1), iv(ierr), iv(ipivot), iv(kalm), p,
     1               v(qtr1), v(h1), v(step1), v, v(w1))
         go to 310
c
c  ***  compute goldfeld-quandt-trotter step (augmented model)  ***
c
 240  if (iv(h) .gt. 0) go to 300
c
c     ***  set h to  d**-1 * ( (j**t)*j + s) ) * d**-1.  ***
c
         h1 = -iv(h)
         iv(h) = h1
         s1 = iv(s)
         if (-1 .ne. iv(kalm)) go to 270
c
c        ***  j is in its original form  ***
c
              do 260 i = 1, p
                   t = one / d(i)
                   do 250 k = 1, i
                        v(h1) = t*(dotprd(n,j(1,i),j(1,k))+v(s1)) / d(k)
                        h1 = h1 + 1
                        s1 = s1 + 1
 250                    continue
 260               continue
              go to 300
c
c  ***  lmstep has applied qrfact to j  ***
c
 270     smh = s1 - h1
         h0 = h1 - 1
         ipiv1 = iv(ipivot)
         t1 = one / d(ipiv1)
         rd0 = iv(rd) - 1
         rdof1 = v(rd0 + 1)
         do 290 i = 1, p
              l = ipiv0 + i
              ipivi = iv(l)
              h1 = h0 + ipivi*(ipivi-1)/2
              l = h1 + ipivi
              m = l + smh
c             ***  v(l) = h(ipivot(i), ipivot(i))  ***
c             ***  v(m) = s(ipivot(i), ipivot(i))  ***
              t = one / d(ipivi)
              rdk = rd0 + i
              e = v(rdk)**2
              if (i .gt. 1) e = e + dotprd(i-1, j(1,i), j(1,i))
              v(l) = (e + v(m)) * t**2
              if (i .eq. 1) go to 290
              l = h1 + ipiv1
              if (ipivi .lt. ipiv1) l = l +
     1                               ((ipiv1-ipivi)*(ipiv1+ipivi-3))/2
              m = l + smh
c             ***  v(l) = h(ipivot(i), ipivot(1))  ***
c             ***  v(m) = s(ipivot(i), ipivot(1))  ***
              v(l) = t * (rdof1 * j(1,i)  +  v(m)) * t1
              if (i .eq. 2) go to 290
              im1 = i - 1
              do 280 k = 2, im1
                   ipk = ipiv0 + k
                   ipivk = iv(ipk)
                   l = h1 + ipivk
                   if (ipivi .lt. ipivk) l = l +
     1                               ((ipivk-ipivi)*(ipivk+ipivi-3))/2
                   m = l + smh
c                  ***  v(l) = h(ipivot(i), ipivot(k))  ***
c                  ***  v(m) = s(ipivot(i), ipivot(k))  ***
                   km1 = k - 1
                   rdk = rd0 + k
                   v(l) = t * (dotprd(km1, j(1,i), j(1,k)) +
     1                            v(rdk)*j(k,i) + v(m)) / d(ipivk)
 280               continue
 290          continue
c
c  ***  compute actual goldfeld-quandt-trotter step  ***
c
 300  h1 = iv(h)
      dig1 = iv(dig)
      lmat1 = iv(lmat)
      call gqtstp(d, v(dig1), v(h1), iv(kagqt), v(lmat1), p, v(step1),
     1            v, v(w1))
c
c
c  ***  compute r(x0 + step)  ***
c
 310  if (iv(irc) .eq. 6) go to 350
      x01 = iv(x0)
      step1 = iv(step)
      call vaxpy(p, x, one, v(step1), v(x01))
      iv(nfcall) = iv(nfcall) + 1
      iv(1) = 1
      iv(toobig) = 0
      go to 999
c
c. . . . . . . . . . . . .  assess candidate step  . . . . . . . . . . .
c
 350  step1 = iv(step)
      lstgst = iv(stlstg)
      x01 = iv(x0)
      call assess(d, iv, p, v(step1), v(lstgst), v, x, v(x01))
c
c  ***  if necessary, switch models and/or restore r  ***
c
      if (iv(switch) .eq. 0) go to 360
         iv(h) = -iabs(iv(h))
         iv(sused) = iv(sused) + 2
         call vcopy(nvsave, v, v(vsave1))
 360  if (iv(restor) .eq. 0) go to 390
         rsave1 = iv(rsave)
         call vcopy(n, r, v(rsave1))
 390  l = iv(irc) - 4
      stpmod = iv(model)
      if (l .gt. 0) go to (410,440,450,450,450,450,450,450,640,570), l
c
c  ***  decide whether to change models  ***
c
      e = v(preduc) - v(fdif)
      sstep = iv(lky)
      s1 = iv(s)
      call slvmul(p, v(sstep), v(s1), v(step1))
      sttsst = half * dotprd(p, v(step1), v(sstep))
      if (iv(model) .eq. 1) sttsst = -sttsst
      if (dabs(e + sttsst) * v(fuzz) .ge. dabs(e)) go to 400
c
c     ***  switch models  ***
c
         iv(model) = 3 - iv(model)
         if (iv(model) .eq. 1) iv(kagqt) = -1
         if (iv(model) .eq. 2 .and. iv(kalm) .gt. 0) iv(kalm) = 0
         if (-2 .lt. l) go to 480
              iv(h) = -iabs(iv(h))
              iv(sused) = iv(sused) + 2
              call vcopy(nvsave, v(vsave1), v)
              go to 420
c
 400  if (-3 .lt. l) go to 480
c
c     ***  recompute step with decreased radius  ***
c
         v(radius) = v(radfac) * v(dstnrm)
         go to 190
c
c  ***  recompute step, saving v values and r if necessary  ***
c
 410  v(radius) = v(radfac) * v(dstnrm)
 420  if (v(f) .ge. v(f0)) go to 190
      rsave1 = iv(rsave)
      call vcopy(n, v(rsave1), r)
      go to 190
c
c  ***  compute step of length v(lmax0) for singular convergence test
c
 440  v(radius) = v(lmax0)
      go to 210
c
c  ***  convergence or false convergence  ***
c
 450  iv(cnvcod) = l
      if (v(f) .ge. v(f0)) go to 700
         if (iv(xirc) .eq. 14) go to 700
              iv(xirc) = 14
c
c. . . . . . . . . . . .  process acceptable step  . . . . . . . . . . .
c
 480  iv(covmat) = 0
c
c  ***  set  lky = (j(x0)**t) * r(x)  ***
c
      lky1 = iv(lky)
      if (iv(kalm) .ge. 0) go to 500
c
c     ***  jacobian has not been modified  ***
c
         do 490 i = 1, p
              v(lky1) = dotprd(n, j(1,i), r)
              lky1 = lky1 + 1
 490          continue
         go to 510
c
c  ***  qrfact has been applied to j.  store copy of r in qtr and  ***
c  ***  apply q to it.                                             ***
c
 500  qtr1 = iv(qtr)
      call vcopy(n, v(qtr1), r)
      call qapply(nn, n, p, j, v(qtr1), iv(ierr))
c
c  ***  multiply top p-vector in qtr by permuted upper triangle    ***
c  ***  stored by qrfact in j and rd.                              ***
c
      rd1 = iv(rd)
      temp1 = iv(stlstg)
      call rptmul(3, iv(ipivot), j, nn, p, v(rd1), v(qtr1), v(lky1),
     1            v(temp1))
c
c  ***  see whether to set v(radfac) by gradient tests  ***
c
 510  if (iv(irc) .ne. 3) go to 560
         step1 = iv(step)
         temp1 = iv(stlstg)
         temp2 = iv(x0)
c
c     ***  set  temp1 = hessian * step  for use in gradient tests  ***
c
         if (stpmod .eq. 2) go to 530
c
c        ***  step computed using gauss-newton model  ***
c        ***  -- qrfact has been applied to j         ***
c
              rd1 = iv(rd)
              call rptmul(2, iv(ipivot), j, nn, p, v(rd1),
     1                    v(step1), v(temp1), v(temp2))
              go to 560
c
c     ***  step computed using augmented model  ***
c
 530     h1 = iv(h)
         k = temp2
         do 540 i = 1, p
              v(k) = d(i) * v(step1)
              k = k + 1
              step1 = step1 + 1
 540          continue
         call slvmul(p, v(temp1), v(h1), v(temp2))
         do 550 i = 1, p
              v(temp1) = d(i) * v(temp1)
              temp1 = temp1 + 1
 550          continue
c
c  ***  save old gradient and compute new one  ***
c
 560  iv(ngcall) = iv(ngcall) + 1
      g1 = iv(g)
      g01 = iv(w)
      call vcopy(p, v(g01), v(g1))
      iv(1) = 2
      go to 999
c
c  ***  initializations -- g0 = g - g0, etc.  ***
c
 570  g01 = iv(w)
      g1 = iv(g)
      call vaxpy(p, v(g01), negone, v(g01), v(g1))
      step1 = iv(step)
      temp1 = iv(stlstg)
      temp2 = iv(x0)
      if (iv(irc) .ne. 3) go to 600
c
c  ***  set v(radfac) by gradient tests  ***
c
c     ***  set  temp1 = d**-1 * (hessian * step  +  (g(x0) - g(x)))  ***
c
         k = temp1
         l = g01
         do 580 i = 1, p
              v(k) = (v(k) - v(l)) / d(i)
              k = k + 1
              l = l + 1
 580          continue
c
c        ***  do gradient tests  ***
c
         if (v2norm(p, v(temp1)) .le. v(dgnorm) * v(tuner4))  go to 590
              if (dotprd(p, v(g1), v(step1))
     1                  .ge. v(gtstep) * v(tuner5))  go to 600
 590               v(radfac) = v(incfac)
c
c  ***  finish computing lky = ((j(x) - j(x0))**t) * r  ***
c
c     ***  currently lky = (j(x0)**t) * r  ***
c
 600  lky1 = iv(lky)
      call vaxpy(p, v(lky1), negone, v(lky1), v(g1))
c
c  ***  determine sizing factor v(size)  ***
c
c     ***  set temp1 = s * step  ***
      s1 = iv(s)
      call slvmul(p, v(temp1), v(s1), v(step1))
c
      t1 = dabs(dotprd(p, v(step1), v(temp1)))
      t = dabs(dotprd(p, v(step1), v(lky1)))
      v(size) = one
      if (t .lt. t1) v(size) = t / t1
c
c  ***  update s  ***
c
      call slupdt(v(s1), v(cosmin), p, v(size), v(step1), v(temp1),
     1            v(temp2), v(g01), v(wscale), v(lky1))
      iv(1) = 2
      go to 150
c
c. . . . . . . . . . . . . .  misc. details  . . . . . . . . . . . . . .
c
c  ***  bad parameters to assess  ***
c
 640  iv(1) = 14
      go to 900
c
c  ***  convergence obtained -- compute covariance matrix if desired ***
c
 700  if (iv(covreq) .eq. 0 .and. iv(covprt) .eq. 0) go to 760
      if (iv(covmat) .ne. 0) go to 760
      if (iv(cnvcod) .ge. 7) go to 760
      iv(mode) = 0
 710  call covclc(i, d, iv, j, n, nn, p, r, v, x)
      go to (720, 720, 740, 750), i
 720  iv(nfcov) = iv(nfcov) + 1
      iv(nfcall) = iv(nfcall) + 1
      iv(restor) = i
      iv(1) = 1
      go to 999
c
 730  if (iv(restor) .eq. 1 .or. iv(toobig) .ne. 0) go to 710
      iv(nfgcal) = iv(nfcall)
 740  iv(ngcov) = iv(ngcov) + 1
      iv(ngcall) = iv(ngcall) + 1
      iv(1) = 2
      go to 999
c
 750  iv(mode) = 0
      if (iv(niter) .eq. 0) iv(mode) = -1
c
 760  iv(1) = iv(cnvcod)
      iv(cnvcod) = 0
c
c  ***  print summary of final iteration and other requested items  ***
c
 900  call itsmry(d, iv, p, v, x)
c
 999  return
c
c  ***  last card of nl2itr follows  ***
      end
