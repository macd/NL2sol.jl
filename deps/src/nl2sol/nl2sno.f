      subroutine nl2sno(n, p, x, calcr, iv, v, uiparm, urparm, ufparm)
c
c  ***  like nl2sol, but without calcj -- minimize nonlinear sum of  ***
c  ***  squares using finite-difference jacobian approximations      ***
c  ***  (nl2sol version 2.2)  ***
c
      integer n, p, iv(1), uiparm(1)
      double precision x(p), v(1), urparm(1)
c     dimension iv(60+p),  v(93 + n*p + 3*n + p*(3*p+33)/2)
      external calcr, ufparm
c
c-----------------------------  discussion  ----------------------------
c
c        the parameters for nl2sno are the same as those for nl2sol
c     (which see), except that calcj is omitted.  instead of calling
c     calcj to obtain the jacobian matrix of r at x, nl2sno computes
c     an approximation to it by finite (forward) differences -- see
c     v(dltfdj) below.  nl2sno uses function values only when comput-
c     the covariance matrix (rather than the functions and gradients
c     that nl2sol may use).  to do so, nl2sno sets iv(covreq) to -1 if
c     iv(covprt) = 1 with iv(covreq) = 0 and to minus its absolute
c     value otherwise.  thus v(delta0) is never referenced and only
c     v(dltfdc) matters -- see nl2sol for a description of v(dltfdc).
c        the number of extra calls on calcr used in computing the jaco-
c     bian approximation are not included in the function evaluation
c     count iv(nfcall) and are not otherwise reported.
c
c v(dltfdj)... v(36) helps choose the step size used when computing the
c             finite-difference jacobian matrix.  for differences in-
c             volving x(i), the step size first tried is
c                       v(dltfdj) * max(abs(x(i)), 1/d(i)),
c             where d is the current scale vector (see ref. 1).  (if
c             this step is too big, i.e., if calcr sets nf to 0, then
c             smaller steps are tried until the step size is shrunk be-
c             low 1000 * machep, where machep is the unit roundoff.
c             default = machep**0.5.
c
c  ***  references  ***
c
c 1.  dennis, j.e., gay, d.m., and welsch, r.e. (1981), an adaptive
c             nonlinear least-squares algorithm, acm trans. math.
c             software, vol. 7, no. 3.
c
c  ***  general  ***
c
c     coded by david m. gay.
c     this subroutine was written in connection with research
c     supported by the national science foundation under grants
c     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and
c     mcs-7906671.
c
c+++++++++++++++++++++++++++  declarations  ++++++++++++++++++++++++++++
c
c  ***  intrinsic functions  ***
c/+
      integer iabs
      double precision dabs, dmax1
c/
c  ***  external functions and subroutines  ***
c
      external dfault, itsmry, nl2itr, rmdcon, vscopy
      double precision rmdcon
c
c dfault... supplies default parameter values.
c itsmry... prints iteration summary and info about initial and final x.
c nl2itr... reverse-communication routine that carries out nl2sol algo-
c             rithm.
c rmdcon... returns machine-dependent constants.
c vscopy... sets all elements of a vector to a scalar.
c
      logical strted
      integer dk, d1, i, j1, j1k, k, nf, rn, r1, dinit
      double precision h, hfac, hlim, negpt5, one, xk, zero
c
c  ***  subscripts for iv and v  ***
c
      integer covprt, covreq, d, dltfdj, dtype, j, nfcall, nfgcal, r,
     1        toobig
c
c/6
c      data hfac/1.d+3/, negpt5/-0.5d+0/, one/1.d+0/, zero/0.d+0/
c/7
      parameter (hfac=1.d+3, negpt5=-0.5d+0, one=1.d+0, zero=0.d+0)
c/
c
c  ***  iv subscript values  ***
c
c/6
c      data covprt/14/, covreq/15/, d/27/, dtype/16/, j/33/,
c     1     nfcall/6/, nfgcal/7/, r/50/, toobig/2/
c/7
      parameter (covprt=14, covreq=15, d=27, dtype=16, j=33,
     1     nfcall=6, nfgcal=7, r=50, toobig=2)
c/
c
c  ***  v subscript values  ***
c
c/6
c      data dltfdj/36/, dinit/38/
c/7
      parameter (dltfdj=36)
      save hlim
c/
      data hlim/0.d+0/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
      d1 = 94 + 2*n + p*(3*p + 31)/2
      iv(d) = d1
      r1 = d1 + p
      iv(r) = r1
      j1 = r1 + n
      iv(j) = j1
      rn = j1 - 1
      if (iv(1) .eq. 0) call dfault(iv, v)
      iv(covreq) = -iabs(iv(covreq))
      if (iv(covprt) .ne. 0 .and. iv(covreq) .eq. 0) iv(covreq) = -1
      strted = .true.
      if (iv(1) .ne. 12) go to 80
         strted = .false.
         iv(nfcall) = 1
         iv(nfgcal) = 1
c        ***  initialize scale vector d to ones for computing
c        ***  initial jacobian.
         if (iv(dtype) .gt. 0) call vscopy(p, v(d1), one)
c       if (v(dinit).gt.zero) call vscopy(p, v(d1), v(dinit))
       if (v(dinit).gt.rmdcon(4)) call vscopy(p, v(d1), v(dinit))
c
 10   nf = iv(nfcall)
      call calcr(n, p, x, nf, v(r1), uiparm, urparm, ufparm)
      if (strted) go to 20
         if (nf .gt. 0) go to 30
              iv(1) = 13
              go to 90
c
 20   if (nf .le. 0) iv(toobig) = 1
      go to 80
c
c  ***  compute finite-difference jacobian  ***
c
 30   j1k = j1
      dk = d1
      do 70 k = 1, p
         xk = x(k)
         h = v(dltfdj) * dmax1(dabs(xk), one/v(dk))
         dk = dk + 1
 40      x(k) = xk + h
         nf = iv(nfgcal)
         call calcr (n, p, x, nf, v(j1k), uiparm, urparm, ufparm)
         if (nf .gt. 0) go to 50
              if (hlim .eq. zero) hlim = hfac * rmdcon(3)
c             ***  hlim = hfac times the unit roundoff  ***
              h = negpt5 * h
              if (dabs(h) .ge. hlim) go to 40
                   iv(1) = 15
                   go to 90
 50      x(k) = xk
         do 60 i = r1, rn
              v(j1k) = (v(j1k) - v(i)) / h
              j1k = j1k + 1
 60           continue
 70      continue
c
      strted = .true.
c
 80   call nl2itr(v(d1), iv, v(j1), n, n, p, v(r1), v, x)
      if (iv(1) - 2) 10, 30, 999
c
 90   call itsmry(v(d1), iv, p, v, x)
c
 999  return
c  ***  last card of nl2sno follows  ***
      end
