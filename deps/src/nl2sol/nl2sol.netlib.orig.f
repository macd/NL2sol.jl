c  algorithm 573
c
c  nl2sol -- an adaptive nonlinear least-squares algorithm
c
c  authors = john e. dennis, jr., david m. gay, and roy e. welsch
c
c  acm transactions on mathematical software, september, 1981.
c
c  this file comes in 9 sections, separated by a comment line having c
c  in column 1 and slashes in columns 2-72.  the first section con-
c  sists of these comments.  sections 2-5 contain single-precision 1966
c  ansi standard fortran source code, and sections 6-9 are double-
c  precision versions of sections 2-5.  comments in sections 4 and 8
c  describe an easy way to modify this code for use with fortran 77.
c  the 9 sections are as follows...
c
c     1. these comments.
c     2. single-prec. short test program.
c     3. single-prec. machine-dependent functions imdcon and rmdcon.
c     4. single-prec. machine-independent nl2sol modules.
c     5. single-prec. long test program.
c     6. double-prec. short test program.
c     7. double-prec. machine-dependent functions imdcon and rmdcon.
c     8. double-prec. machine-independent nl2sol modules.
c     9. double-prec. long test program.
c
c  the short test program (sections 2 and 6) amounts to the example in
c  section 3.2 of the description of toms algorithm 573 with an added
c  call on nl2sno.
c
c  depending on the computer used, it may be necessary to change the
c  data statements in sections 3 and 7 -- see section 3.12 of the
c  description of toms algorithm 573.  (the version of rmdcon in
c  section 3 is set for cdc computers, and that in section 8 is set for
c  ibm 360 and 370 computers.)
c
c  the first three modules in sections 4 and 8 are nl2sol, nl2sno, and
c  nl2itr.  the remaining modules follow in alphabetical order.
c
c  the long test program (sections 5 and 9) runs the tests reported in
c  table ii of the toms paper on nl2sol.  this program produces a
c  one-page summary on unit imdcon(2) and detailed output on unit
c  imdcon(1).  the latter may be suppressed by arranging for imdcon(1)
c  to return 0.
c
c///////////////////////////////////////////////////////////////////////
      end
c     ***  test nl2sol and nl2sno on madsen example  ***
      integer iv(62), uiparm(1)
      real v(147), x(2), urparm(1)
      external madr, madj
      x(1) = 3.0
      x(2) = 1.0
      iv(1) = 0
      call nl2sol(3, 2, x, madr, madj, iv, v, uiparm, urparm, madr)
      iv(1) = 12
      x(1) = 3.0
      x(2) = 1.0
      call nl2sno(3, 2, x, madr, iv, v, uiparm, urparm, madr)
      stop
      end
      subroutine madr(n, p, x, nf, r, uiparm, urparm, ufparm)
      integer n, p, nf, uiparm(1)
      real x(p), r(n), urparm(1)
      external ufparm
      r(1) = x(1)**2 + x(2)**2 + x(1)*x(2)
      r(2) = sin(x(1))
      r(3) = cos(x(2))
      return
      end
      subroutine madj(n, p, x, nf, j, uiparm, urparm, ufparm)
      integer n, p, nf, uiparm(1)
      real x(p), j(n,p), urparm(1)
      external ufparm
      j(1,1) = 2.0*x(1) + x(2)
      j(1,2) = 2.0*x(2) + x(1)
      j(2,1) = cos(x(1))
      j(2,2) = 0.0
      j(3,1) = 0.0
      j(3,2) = -sin(x(2))
      return
      end
      integer function imdcon(k)
c
      integer k
c
c  ***  return integer machine-dependent constants  ***
c
c     ***  k = 1 means return standard output unit number.   ***
c     ***  k = 2 means return alternate output unit number.  ***
c     ***  k = 3 means return  input unit number.            ***
c          (note -- k = 2, 3 are used only by test programs.)
c
      integer mdcon(3)
      data mdcon(1)/6/, mdcon(2)/8/, mdcon(3)/5/
c
      imdcon = mdcon(k)
      return
c  ***  last card of imdcon follows  ***
      end
      real function rmdcon(k)
c
c  ***  return machine dependent constants used by nl2sol  ***
c
c +++  comments below contain data statements for various machines.  +++
c +++  to convert to another machine, place a c in column 1 of the   +++
c +++  data statement line(s) that correspond to the current machine +++
c +++  and remove the c from column 1 of the data statement line(s)  +++
c +++  that correspond to the new machine.                           +++
c
      integer k
c
c  ***  the constant returned depends on k...
c
c  ***        k = 1... smallest pos. eta such that -eta exists.
c  ***        k = 2... square root of 1.001*eta.
c  ***        k = 3... unit roundoff = smallest pos. no. machep such
c  ***                 that 1 + machep .gt. 1 .and. 1 - machep .lt. 1.
c  ***        k = 4... square root of 0.999*machep.
c  ***        k = 5... square root of 0.999*big (see k = 6).
c  ***        k = 6... largest machine no. big such that -big exists.
c
      real big, eta, machep
c/+
      real sqrt
c/
      real one001, pt999
c
      data one001/1.001/, pt999/0.999/
c
c  +++  ibm 360, ibm 370, or xerox  +++
c
c     data big/z7fffffff/, eta/z00100000/, machep/z3c100000/
c
c  +++  data general  +++
c
c     data big/0.7237e+76/, eta/0.5398e-78/, machep/0.9537e-06/
c
c  +++  dec 11  +++
c
c     data big/1.7e+38/, eta/2.9388e-39/, machep/1.1921e-07/
c
c  +++  hp3000  +++
c
c     data big/1.1579e+77/, eta/8.6362e-78/, machep/2.3842e-07/
c
c  +++  honeywell  +++
c
c     data big/o376777000000/, eta/o404400400000/,
c    1     machep/o716400000000/
c
c  +++  dec10  +++
c
c     data big/"377777777777/, eta/"000400000021/,
c    1     machep/"147400000000/
c
c  +++  burroughs  +++
c
c     data big/o0777777777777777/, eta/o1771000000000000/,
c    1     machep/o1301000000000000/
c
c  +++  control data  +++
c
      data big/37754000000000000000b/, eta/00024000000000000000b/,
     1     machep/16414000000000000000b/
c
c  +++  prime  +++
c
c     data big/1.7e+38/, eta/1.47e-39/, machep/2.38419e-7/
c
c  +++  univac  +++
c
c     data big/1.69e+38/, eta/5.9e-39/, machep/1.4901162e-8/
c
c  +++  vax  +++
c
c     data big/1.7e+38/, eta/2.939e-39/, machep/5.9604645e-08/
c
c-------------------------------  body  --------------------------------
c
      go to (10, 20, 30, 40, 50, 60), k
c
 10   rmdcon = eta
      go to 999
c
 20   rmdcon = sqrt(one001*eta)
      go to 999
c
 30   rmdcon = machep
      go to 999
c
 40   rmdcon = sqrt(pt999*machep)
      go to 999
c
 50   rmdcon = sqrt(pt999*big)
      go to 999
c
 60   rmdcon = big
c
 999  return
c  ***  last card of rmdcon follows  ***
      end
      subroutine nl2sol(n, p, x, calcr, calcj, iv, v, uiparm, urparm,
     1                  ufparm)
c
c  ***  minimize nonlinear sum of squares using analytic jacobian  ***
c  ***  (nl2sol version 2.2)  ***
c
      integer n, p, iv(1), uiparm(1)
      real x(p), v(1), urparm(1)
c     dimension iv(60+p),  v(93 + n*p + 3*n + p*(3*p+33)/2)
c     dimension uiparm(*), urparm(*)
      external calcr, calcj, ufparm
c
c  ***  purpose  ***
c
c        given a p-vector x of parameters, calcr computes an n-vector
c     r = r(x) of residuals corresponding to x.  (r(x) probably arises
c     from a nonlinear model involving p parameters and n observations.)
c     this routine interacts with nl2itr to seek a parameter vector x
c     that minimizes the sum of the squares of (the components of) r(x),
c     i.e., that minimizes the sum-of-squares function
c     f(x) = (r(x)**t) * r(x) / 2.  r(x) is assumed to be a twice con-
c     tinuously differentiable function of x.
c
c--------------------------  parameter usage  --------------------------
c
c n........ (input) the number of observations, i.e., the number of
c                  components in r(x).  n must be .ge. p.
c p........ (input) the number of parameters (components in x).  p must
c                  be positive.
c x........ (input/output).  on input, x is an initial guess at the
c                  desired parameter estimate.  on output, x contains
c                  the best parameter estimate found.
c calcr.... (input) a subroutine which, given x, computes r(x).  calcr
c                  must be declared external in the calling program.
c                  it is invoked by
c                       call calcr(n,p,x,nf,r,uiparm,urparm,ufparm)
c                  when calcr is called, nf is the invocation count
c                  for calcr.  it is included for possible use with
c                  calcj.  if x is out of bounds (e.g. if it would
c                  cause overflow in computing r(x)), then calcr should
c                  set nf to 0.  this will cause a shorter step to be
c                  attempted.  the other parameters are as described
c                  above and below.  calcr should not change n, p, or x.
c calcj.... (input) a subroutine which, given x, computes the jacobian
c                  matrix j of r at x, i.e., the n by p matrix whose
c                  (i,k) entry is the partial derivative of the i-th
c                  component of r with respect to x(k).  calcj must be
c                  declared external in the calling program.  it is
c                  invoked by
c                       call calcj(n,p,x,nf,j,uiparm,urparm,ufparm)
c                  nf is the invocation count for calcr at the time
c                  r(x) was evaluated.  the x passed to calcj is
c                  usually the one passed to calcr on either its most
c                  recent invocation or the one prior to it.  if calcr
c                  saves intermediate results for use by calcj, then it
c                  is possible to tell from nf whether they are valid
c                  for the current x (or which copy is valid if two
c                  copies are kept).  if j cannot be computed at x,
c                  then calcj should set nf to 0.  in this case, nl2sol
c                  will return with iv(1) = 15.  the other parameters
c                  to calcj are as described above and below.  calcj
c                  should not change n, p, or x.
c iv....... (input/output) an integer value array of length at least
c                  60 + p that helps control the nl2sol algorithm and
c                  that is used to store various intermediate quanti-
c                  ties.  of particular interest are the initialization/
c                  return code iv(1) and the entries in iv that control
c                  printing and limit the number of iterations and func-
c                  tion evaluations.  see the section on iv input
c                  values below.
c v........ (input/output) a floating-point value array of length at
c                  least 93 + n*p + 3*n + p*(3*p+33)/2 that helps con-
c                  trol the nl2sol algorithm and that is used to store
c                  various intermediate quantities.  of particular in-
c                  terest are the entries in v that limit the length of
c                  the first step attempted (lmax0), specify conver-
c                  gence tolerances (afctol, rfctol, xctol, xftol),
c                  and help choose the step size used in computing the
c                  covariance matrix (delta0).  see the section on
c                  (selected) v input values below.
c uiparm... (input) user integer parameter array passed without change
c                  to calcr and calcj.
c urparm... (input) user floating-point parameter array passed without
c                  change to calcr and calcj.
c ufparm... (input) user external subroutine or function passed without
c                  change to calcr and calcj.
c
c  ***  iv input values (from subroutine dfault)  ***
c
c iv(1)...  on input, iv(1) should have a value between 0 and 12......
c             0 and 12 mean this is a fresh start.  0 means that
c             dfault(iv, v) is to be called to provide all default
c             values to iv and v.  12 (the value that dfault assigns to
c             iv(1)) means the caller has already called dfault(iv, v)
c             and has possibly changed some iv and/or v entries to non-
c             default values.  default = 12.
c iv(covprt)... iv(14) = 1 means print a covariance matrix at the solu-
c             tion.  (this matrix is computed just before a return with
c             iv(1) = 3, 4, 5, 6.)
c             iv(covprt) = 0 means skip this printing.  default = 1.
c iv(covreq)... iv(15) = nonzero means compute a covariance matrix
c             just before a return with iv(1) = 3, 4, 5, 6.  in
c             this case, an approximate covariance matrix is obtained
c             in one of several ways.  let k = abs(iv(covreq)) and let
c             scale = 2*f(x)/max(1,n-p),  where 2*f(x) is the residual
c             sum of squares.  if k = 1 or 2, then a finite-difference
c             hessian approximation h is obtained.  if h is positive
c             definite (or, for k = 3, if the jacobian matrix j at x
c             is nonsingular), then one of the following is computed...
c                  k = 1....  scale * h**-1 * (j**t * j) * h**-1.
c                  k = 2....  scale * h**-1.
c                  k = 3....  scale * (j**t * j)**-1.
c             (j**t is the transpose of j, while **-1 means inverse.)
c             if iv(covreq) is positive, then both function and grad-
c             ient values (calls on calcr and calcj) are used in com-
c             puting h (with step sizes determined using v(delta0) --
c             see below), while if iv(covreq) is negative, then only
c             function values (calls on calcr) are used (with step
c             sizes determined using v(dltfdc) -- see below).  if
c             iv(covreq) = 0, then no attempt is made to compute a co-
c             variance matrix (unless iv(covprt) = 1, in which case
c             iv(covreq) = 1 is assumed).  see iv(covmat) below.
c             default = 1.
c iv(dtype).... iv(16) tells how the scale vector d (see ref. 1) should
c             be chosen.  iv(dtype) .ge. 1 means choose d as described
c             below with v(dfac).  iv(dtype) .le. 0 means the caller
c             has chosen d and has stored it in v starting at
c             v(94 + 2*n + p*(3*p + 31)/2).  default = 1.
c iv(inits).... iv(25) tells how the s matrix (see ref. 1) should be
c             initialized.  0 means initialize s to 0 (and start with
c             the gauss-newton model).  1 and 2 mean that the caller
c             has stored the lower triangle of the initial s rowwise in
c             v starting at v(87+2*p).  iv(inits) = 1 means start with
c             the gauss-newton model, while iv(inits) = 2 means start
c             with the augmented model (see ref. 1).  default = 0.
c iv(mxfcal)... iv(17) gives the maximum number of function evaluations
c             (calls on calcr, excluding those used to compute the co-
c             variance matrix) allowed.  if this number does not suf-
c             fice, then nl2sol returns with iv(1) = 9.  default = 200.
c iv(mxiter)... iv(18) gives the maximum number of iterations allowed.
c             it also indirectly limits the number of gradient evalua-
c             tions (calls on calcj, excluding those used to compute
c             the covariance matrix) to iv(mxiter) + 1.  if iv(mxiter)
c             iterations do not suffice, then nl2sol returns with
c             iv(1) = 10.  default = 150.
c iv(outlev)... iv(19) controls the number and length of iteration sum-
c             mary lines printed (by itsmry).  iv(outlev) = 0 means do
c             not print any summary lines.  otherwise, print a summary
c             line after each abs(iv(outlev)) iterations.  if iv(outlev)
c             is positive, then summary lines of length 117 (plus carri-
c             age control) are printed, including the following...  the
c             iteration and function evaluation counts, current func-
c             tion value (v(f) = half the sum of squares), relative
c             difference in function values achieved by the latest step
c             (i.e., reldf = (f0-v(f))/f0, where f0 is the function
c             value from the previous iteration), the relative function
c             reduction predicted for the step just taken (i.e.,
c             preldf = v(preduc) / f0, where v(preduc) is described
c             below), the scaled relative change in x (see v(reldx)
c             below), the models used in the current iteration (g =
c             gauss-newton, s=augmented), the marquardt parameter
c             stppar used in computing the last step, the sizing factor
c             used in updating s, the 2-norm of the scale vector d
c             times the step just taken (see ref. 1), and npreldf, i.e.,
c             v(nreduc)/f0, where v(nreduc) is described below -- if
c             npreldf is positive, then it is the relative function
c             reduction predicted for a newton step (one with
c             stppar = 0).  if npreldf is zero, either the gradient
c             vanishes (as does preldf) or else the augmented model
c             is being used and its hessian is indefinite (with preldf
c             positive).  if npreldf is negative, then it is the nega-
c             of the relative function reduction predicted for a step
c             computed with step bound v(lmax0) for use in testing for
c             singular convergence.
c                  if iv(outlev) is negative, then lines of maximum
c             length 79 (or 55 is iv(covprt) = 0) are printed, includ-
c             ing only the first 6 items listed above (through reldx).
c             default = 1.
c iv(parprt)... iv(20) = 1 means print any nondefault v values on a
c             fresh start or any changed v values on a restart.
c             iv(parprt) = 0 means skip this printing.  default = 1.
c iv(prunit)... iv(21) is the output unit number on which all printing
c             is done.  iv(prunit) = 0 means suppress all printing.
c             (setting iv(prunit) to 0 is the only way to suppress the
c             one-line termination reason message printed by itsmry.)
c             default = standard output unit (unit 6 on most systems).
c iv(solprt)... iv(22) = 1 means print out the value of x returned (as
c             well as the corresponding gradient and scale vector d).
c             iv(solprt) = 0 means skip this printing.  default = 1.
c iv(statpr)... iv(23) = 1 means print summary statistics upon return-
c             ing.  these consist of the function value (half the sum
c             of squares) at x, v(reldx) (see below), the number of
c             function and gradient evaluations (calls on calcr and
c             calcj respectively, excluding any calls used to compute
c             the covariance), the relative function reductions predict-
c             ed for the last step taken and for a newton step (or per-
c             haps a step bounded by v(lmax0) -- see the descriptions
c             of preldf and npreldf under iv(outlev) above), and (if an
c             attempt was made to compute the covariance) the number of
c             calls on calcr and calcj used in trying to compute the
c             covariance.  iv(statpr) = 0 means skip this printing.
c             default = 1.
c iv(x0prt).... iv(24) = 1 means print the initial x and scale vector d
c             (on a fresh start only).  iv(x0prt) = 0 means skip this
c             printing.  default = 1.
c
c  ***  (selected) iv output values  ***
c
c iv(1)........ on output, iv(1) is a return code....
c             3 = x-convergence.  the scaled relative difference be-
c                  tween the current parameter vector x and a locally
c                  optimal parameter vector is very likely at most
c                  v(xctol).
c             4 = relative function convergence.  the relative differ-
c                  ence between the current function value and its lo-
c                  cally optimal value is very likely at most v(rfctol).
c             5 = both x- and relative function convergence (i.e., the
c                  conditions for iv(1) = 3 and iv(1) = 4 both hold).
c             6 = absolute function convergence.  the current function
c                  value is at most v(afctol) in absolute value.
c             7 = singular convergence.  the hessian near the current
c                  iterate appears to be singular or nearly so, and a
c                  step of length at most v(lmax0) is unlikely to yield
c                  a relative function decrease of more than v(rfctol).
c             8 = false convergence.  the iterates appear to be converg-
c                  ing to a noncritical point.  this may mean that the
c                  convergence tolerances (v(afctol), v(rfctol),
c                  v(xctol)) are too small for the accuracy to which
c                  the function and gradient are being computed, that
c                  there is an error in computing the gradient, or that
c                  the function or gradient is discontinuous near x.
c             9 = function evaluation limit reached without other con-
c                  vergence (see iv(mxfcal)).
c            10 = iteration limit reached without other convergence
c                  (see iv(mxiter)).
c            11 = stopx returned .true. (external interrupt).  see the
c                  usage notes below.
c            13 = f(x) cannot be computed at the initial x.
c            14 = bad parameters passed to assess (which should not
c                  occur).
c            15 = the jacobian could not be computed at x (see calcj
c                  above).
c            16 = n or p (or parameter nn to nl2itr) out of range --
c                  p .le. 0 or n .lt. p or nn .lt. n.
c            17 = restart attempted with n or p (or par. nn to nl2itr)
c                  changed.
c            18 = iv(inits) is out of range.
c            19...45 = v(iv(1)) is out of range.
c            50 = iv(1) was out of range.
c            87...(86+p) = jtol(iv(1)-86) (i.e., v(iv(1)) is not
c                  positive (see v(dfac) below).
c iv(covmat)... iv(26) tells whether a covariance matrix was computed.
c             if (iv(covmat) is positive, then the lower triangle of
c             the covariance matrix is stored rowwise in v starting at
c             v(iv(covmat)).  if iv(covmat) = 0, then no attempt was
c             made to compute the covariance.  if iv(covmat) = -1,
c             then the finite-difference hessian was indefinite.  and
c             and if iv(covmat) = -2, then a successful finite-differ-
c             encing step could not be found for some component of x
c             (i.e., calcr set nf to 0 for each of two trial steps).
c             note that iv(covmat) is reset to 0 after each successful
c             step, so if such a step is taken after a restart, then
c             the covariance matrix will be recomputed.
c iv(d)........ iv(27) is the starting subscript in v of the current
c             scale vector d.
c iv(g)........ iv(28) is the starting subscript in v of the current
c             least-squares gradient vector (j**t)*r.
c iv(nfcall)... iv(6) is the number of calls so far made on calcr (i.e.,
c             function evaluations, including those used in computing
c             the covariance).
c iv(nfcov).... iv(40) is the number of calls made on calcr when
c             trying to compute covariance matrices.
c iv(ngcall)... iv(30) is the number of gradient evaluations (calls on
c             calcj) so far done (including those used for computing
c             the covariance).
c iv(ngcov).... iv(41) is the number of calls made on calcj when
c             trying to compute covariance matrices.
c iv(niter).... iv(31) is the number of iterations performed.
c iv(r)........ iv(50) is the starting subscript in v of the residual
c             vector r corresponding to x.
c
c  ***  (selected) v input values (from subroutine dfault)  ***
c
c v(afctol)... v(31) is the absolute function convergence tolerance.
c             if nl2sol finds a point where the function value (half
c             the sum of squares) is less than v(afctol), and if nl2sol
c             does not return with iv(1) = 3, 4, or 5, then it returns
c             with iv(1) = 6.  default = max(10**-20, machep**2), where
c             machep is the unit roundoff.
c v(delta0)... v(44) is a factor used in choosing the finite-difference
c             step size used in computing the covariance matrix when
c             iv(covreq) = 1 or 2.  for component i, step size
c                  v(delta0) * max(abs(x(i)), 1/d(i)) * sign(x(i))
c             is used, where d is the current scale vector (see ref. 1).
c             (if this step results in calcr setting nf to 0, then -0.5
c             times this step is also tried.)  default = machep**0.5,
c             where machep is the unit roundoff.
c v(dfac)..... v(41) and the d0 and jtol arrays (see v(d0init) and
c             v(jtinit)) are used in updating the scale vector d when
c             iv(dtype) .gt. 0.  (d is initialized according to
c             v(dinit).)  let d1(i) =
c               max(sqrt(jcnorm(i)**2 + max(s(i,i),0)), v(dfac)*d(i)),
c             where jcnorm(i) is the 2-norm of the i-th column of the
c             current jacobian matrix and s is the s matrix of ref. 1.
c             if iv(dtype) = 1, then d(i) is set to d1(i) unless
c             d1(i) .lt. jtol(i), in which case d(i) is set to
c                                max(d0(i), jtol(i)).
c             if iv(dtype) .ge. 2, then d is updated during the first
c             iteration as for iv(dtype) = 1 (after any initialization
c             due to v(dinit)) and is left unchanged thereafter.
c             default = 0.6.
c v(dinit).... v(38), if nonnegative, is the value to which the scale
c             vector d is initialized.  default = 0.
c v(dltfdc)... v(40) helps choose the step size used when computing the
c             covariance matrix when iv(covreq) = -1 or -2.  for
c             differences involving x(i), the step size first tried is
c                       v(dltfdc) * max(abs(x(i)), 1/d(i)),
c             where d is the current scale vector (see ref. 1).  (if
c             this step is too big the first time it is tried, i.e., if
c             calcr sets nf to 0, then -0.5 times this step is also
c             tried.)  default = machep**(1/3), where machep is the
c             unit roundoff.
c v(d0init)... v(37), if positive, is the value to which all components
c             of the d0 vector (see v(dfac)) are initialized.  if
c             v(dfac) = 0, then it is assumed that the caller has
c             stored d0 in v starting at v(p+87).  default = 1.0.
c v(jtinit)... v(39), if positive, is the value to which all components
c             of the jtol array (see v(dfac)) are initialized.  if
c             v(jtinit) = 0, then it is assumed that the caller has
c             stored jtol in v starting at v(87).  default = 10**-6.
c v(lmax0).... v(35) gives the maximum 2-norm allowed for d times the
c             very first step that nl2sol attempts.  it is also used
c             in testing for singular convergence -- if the function
c             reduction predicted for a step of length bounded by
c             v(lmax0) is at most v(rfctol) * abs(f0), where  f0  is
c             the function value at the start of the current iteration,
c             and if nl2sol does not return with iv(1) = 3, 4, 5, or 6,
c             then it returns with iv(1) = 7.    default = 100.
c v(rfctol)... v(32) is the relative function convergence tolerance.
c             if the current model predicts a maximum possible function
c             reduction (see v(nreduc)) of at most v(rfctol)*abs(f0) at
c             the start of the current iteration, where  f0  is the
c             then current function value, and if the last step attempt-
c             ed achieved no more than twice the predicted function
c             decrease, then nl2sol returns with iv(1) = 4 (or 5).
c             default = max(10**-10, machep**(2/3)), where machep is
c             the unit roundoff.
c v(tuner1)... v(26) helps decide when to check for false convergence
c             and to consider switching models.  this is done if the
c             actual function decrease from the current step is no more
c             than v(tuner1) times its predicted value.  default = 0.1.
c v(xctol).... v(33) is the x-convergence tolerance.  if a newton step
c             (see v(nreduc)) is tried that has v(reldx) .le. v(xctol)
c             and if this step yields at most twice the predicted func-
c             tion decrease, then nl2sol returns with iv(1) = 3 (or 5).
c             (see the description of v(reldx) below.)
c             default = machep**0.5, where machep is the unit roundoff.
c v(xftol).... v(34) is the false convergence tolerance.  if a step is
c             tried that gives no more than v(tuner1) times the predict-
c             ed function decrease and that has v(reldx) .le. v(xftol),
c             and if nl2sol does not return with iv(1) = 3, 4, 5, 6, or
c             7, then it returns with iv(1) = 8.  (see the description
c             of v(reldx) below.)  default = 100*machep, where
c             machep is the unit roundoff.
c v(*)........ dfault supplies to v a number of tuning constants, with
c             which it should ordinarily be unnecessary to tinker.  see
c             version 2.2 of the nl2sol usage summary (which is an
c             appendix to ref. 1).
c
c  ***  (selected) v output values  ***
c
c v(dgnorm)... v(1) is the 2-norm of (d**-1)*g, where g is the most re-
c             cently computed gradient and d is the corresponding scale
c             vector.
c v(dstnrm)... v(2) is the 2-norm of d*step, where step is the most re-
c             cently computed step and d is the current scale vector.
c v(f)........ v(10) is the current function value (half the sum of
c             squares).
c v(f0)....... v(13) is the function value at the start of the current
c             iteration.
c v(nreduc)... v(6), if positive, is the maximum function reduction
c             possible according to the current model, i.e., the func-
c             tion reduction predicted for a newton step (i.e.,
c             step = -h**-1 * g,  where  g = (j**t) * r  is the current
c             gradient and h is the current hessian approximation --
c             h = (j**t)*j  for the gauss-newton model and
c             h = (j**t)*j + s  for the augmented model).
c                  v(nreduc) = zero means h is not positive definite.
c                  if v(nreduc) is negative, then it is the negative of
c             the function reduction predicted for a step computed with
c             a step bound of v(lmax0) for use in testing for singular
c             convergence.
c v(preduc)... v(7) is the function reduction predicted (by the current
c             quadratic model) for the current step.  this (divided by
c             v(f0)) is used in testing for relative function
c             convergence.
c v(reldx).... v(17) is the scaled relative change in x caused by the
c             current step, computed as
c                  max(abs(d(i)*(x(i)-x0(i)), 1 .le. i .le. p) /
c                     max(d(i)*(abs(x(i))+abs(x0(i))), 1 .le. i .le. p),
c             where x = x0 + step.
c
c-------------------------------  notes  -------------------------------
c
c  ***  algorithm notes  ***
c
c        see ref. 1 for a description of the algorithm used.
c        on problems which are naturally well scaled, better perform-
c     ance may be obtained by setting v(d0init) = 1.0 and iv(dtype) = 0,
c     which will cause the scale vector d to be set to all ones.
c
c  ***  usage notes  ***
c
c        after a return with iv(1) .le. 11, it is possible to restart,
c     i.e., to change some of the iv and v input values described above
c     and continue the algorithm from the point where it was interrupt-
c     ed.  iv(1) should not be changed, nor should any entries of iv
c     and v other than the input values (those supplied by dfault).
c        those who do not wish to write a calcj which computes the ja-
c     cobian matrix analytically should call nl2sno rather than nl2sol.
c     nl2sno uses finite differences to compute an approximate jacobian.
c        those who would prefer to provide r and j (the residual and
c     jacobian) by reverse communication rather than by writing subrou-
c     tines calcr and calcj may call on nl2itr directly.  see the com-
c     ments at the beginning of nl2itr.
c        those who use nl2sol interactively may wish to supply their
c     own stopx function, which should return .true. if the break key
c     has been pressed since stopx was last invoked.  this makes it pos-
c     sible to externally interrupt nl2sol (which will return with
c     iv(1) = 11 if stopx returns .true.).
c        storage for j is allocated at the end of v.  thus the caller
c     may make v longer than specified above and may allow calcj to use
c     elements of j beyond the first n*p as scratch storage.
c
c  ***  portability notes  ***
c
c        the nl2sol distribution tape contains both single- and double-
c     precision versions of the nl2sol source code, so it should be un-
c     necessary to change precisions.
c        only the functions imdcon and rmdcon contain machine-dependent
c     constants.  to change from one machine to another, it should
c     suffice to change the (few) relevant lines in these functions.
c        intrinsic functions are explicitly declared.  on certain com-
c     puters (e.g. univac), it may be necessary to comment out these
c     declarations.  so that this may be done automatically by a simple
c     program, such declarations are preceded by a comment having c/+
c     in columns 1-3 and blanks in columns 4-72 and are followed by
c     a comment having c/ in columns 1 and 2 and blanks in columns 3-72.
c        the nl2sol source code is expressed in 1966 ansi standard
c     fortran.  it may be converted to fortran 77 by
c     commenting out all lines that fall between a line having c/6 in
c     columns 1-3 and a line having c/7 in columns 1-3 and by removing
c     (i.e., replacing by a blank) the c in column 1 of the lines that
c     follow the c/7 line and preceed a line having c/ in columns 1-2
c     and blanks in columns 3-72.  these changes convert some data
c     statements into parameter statements, convert some variables from
c     real to character*4, and make the data statements that initialize
c     these variables use character strings delimited by primes instead
c     of hollerith constants.  (such variables and data statements
c     appear only in modules itsmry and parchk.  parameter statements
c     appear nearly everywhere.)
c
c  ***  references  ***
c
c 1.  dennis, j.e., gay, d.m., and welsch, r.e. (1981), an adaptive
c             nonlinear least-squares algorithm, acm trans. math.
c             software, vol. 7, no. 3.
c
c
c  ***  general  ***
c
c     coded by david m. gay (winter 1979 - winter 1980).
c     this subroutine was written in connection with research
c     supported by the national science foundation under grants
c     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and
c     mcs-7906671.
c
c----------------------------  declarations  ---------------------------
c
      external itsmry, nl2itr
c itsmry... prints iteration summary and info about initial and final x.
c nl2itr... reverse-communication routine that carries out nl2sol algo-
c             rithm.
c
      logical strted
      integer d1, j1, nf, r1
c
c  ***  subscripts for iv and v  ***
c
      integer d, j, nfcall, nfgcal, r, toobig
c
c  ***  iv subscript values  ***
c
c/6
      data nfcall/6/, nfgcal/7/, toobig/2/
c/7
c     parameter (nfcall=6, nfgcal=7, toobig=2)
c/
c
c  ***  v subscript values  ***
c
c/6
      data d/27/, j/33/, r/50/
c/7
c     parameter (d=27, j=33, r=50)
c/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
      d1 = 94 + 2*n + p*(3*p + 31)/2
      iv(d) = d1
      r1 = d1 + p
      iv(r) = r1
      j1 = r1 + n
      iv(j) = j1
      strted = .true.
      if (iv(1) .ne. 0 .and. iv(1) .ne. 12) go to 40
         strted = .false.
         iv(nfcall) = 1
         iv(nfgcal) = 1
c
 10   nf = iv(nfcall)
      call calcr(n, p, x, nf, v(r1), uiparm, urparm, ufparm)
      if (strted) go to 20
         if (nf .gt. 0) go to 30
              iv(1) = 13
              go to 60
c
 20   if (nf .le. 0) iv(toobig) = 1
      go to 40
c
 30   call calcj(n, p, x, iv(nfgcal), v(j1), uiparm, urparm, ufparm)
      if (iv(nfgcal) .eq. 0) go to 50
      strted = .true.
c
 40   call nl2itr(v(d1), iv, v(j1), n, n, p, v(r1), v, x)
      if (iv(1) - 2) 10, 30, 999
c
 50   iv(1) = 15
 60   call itsmry(v(d1), iv, p, v, x)
c
 999  return
c  ***  last card of nl2sol follows  ***
      end
      subroutine nl2sno(n, p, x, calcr, iv, v, uiparm, urparm, ufparm)
c
c  ***  like nl2sol, but without calcj -- minimize nonlinear sum of  ***
c  ***  squares using finite-difference jacobian approximations      ***
c  ***  (nl2sol version 2.2)  ***
c
      integer n, p, iv(1), uiparm(1)
      real x(p), v(1), urparm(1)
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
      real abs, amax1
c/
c  ***  external functions and subroutines  ***
c
      external dfault, itsmry, nl2itr, rmdcon, vscopy
      real rmdcon
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
      real h, hfac, hlim, negpt5, one, xk, zero
c
c  ***  subscripts for iv and v  ***
c
      integer covprt, covreq, d, dltfdj, dtype, j, nfcall, nfgcal, r,
     1        toobig
c
c/6
      data hfac/1.e+3/, negpt5/-0.5e+0/, one/1.e+0/, zero/0.e+0/
c/7
c     parameter (hfac=1.d+3, negpt5=-0.5d+0, one=1.d+0, zero=0.d+0)
c/
c
c  ***  iv subscript values  ***
c
c/6
      data covprt/14/, covreq/15/, d/27/, dtype/16/, j/33/,
     1     nfcall/6/, nfgcal/7/, r/50/, toobig/2/
c/7
c     parameter (covprt=14, covreq=15, d=27, dtype=16, j=33,
c    1     nfcall=6, nfgcal=7, r=50, toobig=2)
c/
c
c  ***  v subscript values  ***
c
c/6
      data dltfdj/36/, dinit/38/
c/7
c     parameter (dltfdj=36)
c     save hlim
c/
      data hlim/0.e+0/
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
       if (v(dinit).gt.zero) call vscopy(p, v(d1), v(dinit))
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
         h = v(dltfdj) * amax1(abs(xk), one/v(dk))
         dk = dk + 1
 40      x(k) = xk + h
         nf = iv(nfgcal)
         call calcr (n, p, x, nf, v(j1k), uiparm, urparm, ufparm)
         if (nf .gt. 0) go to 50
              if (hlim .eq. zero) hlim = hfac * rmdcon(3)
c             ***  hlim = hfac times the unit roundoff  ***
              h = negpt5 * h
              if (abs(h) .ge. hlim) go to 40
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
      subroutine nl2itr (d, iv, j, n, nn, p, r, v, x)
c
c  ***  carry out nl2sol (nonlinear least-squares) iterations  ***
c  ***  (nl2sol version 2.2)  ***
c
c  ***  parameter declarations  ***
c
      integer iv(1), n, nn, p
      real d(p), j(nn,p), r(n), v(1), x(p)
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
      real e, rdof1, sttsst, t, t1
c
c     ***  constants  ***
c
      real half, negone, one, zero
c
c  ***  intrinsic functions  ***
c/+
      integer iabs
      real abs
c/
c  ***  external functions and subroutines  ***
c
      external assess, covclc, dotprd, dupdat, gqtstp, itsmry, lmstep,
     1         parchk, qapply, qrfact, rptmul, slupdt, slvmul, stopx,
     2         vaxpy, vcopy, vscopy, v2norm
      logical stopx
      real dotprd, v2norm
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
      data cnvcod/34/, covmat/26/, covprt/14/,
     1     covreq/15/, dig/43/, dtype/16/, g/28/, h/44/,
     2     ierr/32/, inits/25/, ipivot/61/, ipiv0/60/,
     3     irc/3/, kagqt/35/, kalm/36/, lky/37/, lmat/58/,
     4     mode/38/, model/5/, mxfcal/17/, mxiter/18/,
     5     nfcall/6/, nfgcal/7/, nfcov/40/, ngcov/41/,
     6     ngcall/30/, niter/31/, qtr/49/,
     7     radinc/8/, rd/51/, restor/9/, rsave/52/, s/53/,
     8     step/55/, stglim/11/, stlstg/56/, sused/57/,
     9     switch/12/, toobig/2/, w/59/, xirc/13/, x0/60/
c/7
c     parameter (cnvcod=34, covmat=26, covprt=14,
c    1     covreq=15, dig=43, dtype=16, g=28, h=44,
c    2     ierr=32, inits=25, ipivot=61, ipiv0=60,
c    3     irc=3, kagqt=35, kalm=36, lky=37, lmat=58,
c    4     mode=38, model=5, mxfcal=17, mxiter=18,
c    5     nfcall=6, nfgcal=7, nfcov=40, ngcov=41,
c    6     ngcall=30, niter=31, qtr=49,
c    7     radinc=8, rd=51, restor=9, rsave=52, s=53,
c    8     step=55, stglim=11, stlstg=56, sused=57,
c    9     switch=12, toobig=2, w=59, xirc=13, x0=60)
c/
c
c  ***  v subscript values  ***
c
c/6
      data cosmin/43/, dgnorm/1/, dinit/38/, dstnrm/2/,
     1     d0init/37/, f/10/, fdif/11/, fuzz/45/,
     2     f0/13/, gtstep/4/, incfac/23/,
     3     jtinit/39/, jtol1/87/, lmax0/35/,
     4     nvsave/9/, phmxfc/21/, preduc/7/,
     5     radfac/16/, radius/8/, rad0/9/, rlimit/42/,
     6     size/47/, stppar/5/, tuner4/29/, tuner5/30/,
     7     vsave1/78/, wscale/48/
c/7
c     parameter (cosmin=43, dgnorm=1, dinit=38, dstnrm=2,
c    1     d0init=37, f=10, fdif=11, fuzz=45,
c    2     f0=13, gtstep=4, incfac=23,
c    3     jtinit=39, jtol1=87, lmax0=35,
c    4     nvsave=9, phmxfc=21, preduc=7,
c    5     radfac=16, radius=8, rad0=9, rlimit=42,
c    6     size=47, stppar=5, tuner4=29, tuner5=30,
c    7     vsave1=78, wscale=48)
c/
c
c
c/6
      data half/0.5e+0/, negone/-1.e+0/, one/1.e+0/, zero/0.e+0/
c/7
c     parameter (half=0.5d+0, negone=-1.d+0, one=1.d+0, zero=0.d+0)
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
      if (abs(e + sttsst) * v(fuzz) .ge. abs(e)) go to 400
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
      t1 = abs(dotprd(p, v(step1), v(temp1)))
      t = abs(dotprd(p, v(step1), v(lky1)))
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
      subroutine assess (d, iv, p, step, stlstg, v, x, x0)
c
c  ***  assess candidate step (nl2sol version 2.2)  ***
c
      integer p, iv(13)
      real d(p), step(p), stlstg(p), v(35), x(p), x0(p)
c
c  ***  purpose  ***
c
c        this subroutine is called by an unconstrained minimization
c     routine to assess the next candidate step.  it may recommend one
c     of several courses of action, such as accepting the step, recom-
c     puting it using the same or a new quadratic model, or halting due
c     to convergence or false convergence.  see the return code listing
c     below.
c
c--------------------------  parameter usage  --------------------------
c
c     iv (i/o) integer parameter and scratch vector -- see description
c             below of iv values referenced.
c      d (in)  scale vector used in computing v(reldx) -- see below.
c      p (in)  number of parameters being optimized.
c   step (i/o) on input, step is the step to be assessed.  it is un-
c             changed on output unless a previous step achieved a
c             better objective function reduction, in which case stlstg
c             will have been copied to step.
c stlstg (i/o) when assess recommends recomputing step even though the
c             current (or a previous) step yields an objective func-
c             tion decrease, it saves in stlstg the step that gave the
c             best function reduction seen so far (in the current itera-
c             tion).  if the recomputed step yields a larger function
c             value, then step is restored from stlstg and
c             x = x0 + step is recomputed.
c      v (i/o) real parameter and scratch vector -- see description
c             below of v values referenced.
c      x (i/o) on input, x = x0 + step is the point at which the objec-
c             tive function has just been evaluated.  if an earlier
c             step yielded a bigger function decrease, then x is
c             restored to the corresponding earlier value.  otherwise,
c             if the current step does not give any function decrease,
c             then x is restored to x0.
c     x0 (in)  initial objective function parameter vector (at the
c             start of the current iteration).
c
c  ***  iv values referenced  ***
c
c    iv(irc) (i/o) on input for the first step tried in a new iteration,
c             iv(irc) should be set to 3 or 4 (the value to which it is
c             set when step is definitely to be accepted).  on input
c             after step has been recomputed, iv(irc) should be
c             unchanged since the previous return of assess.
c                on output, iv(irc) is a return code having one of the
c             following values...
c                  1 = switch models or try smaller step.
c                  2 = switch models or accept step.
c                  3 = accept step and determine v(radfac) by gradient
c                       tests.
c                  4 = accept step, v(radfac) has been determined.
c                  5 = recompute step (using the same model).
c                  6 = recompute step with radius = v(lmax0) but do not
c                       evaulate the objective function.
c                  7 = x-convergence (see v(xctol)).
c                  8 = relative function convergence (see v(rfctol)).
c                  9 = both x- and relative function convergence.
c                 10 = absolute function convergence (see v(afctol)).
c                 11 = singular convergence (see v(lmax0)).
c                 12 = false convergence (see v(xftol)).
c                 13 = iv(irc) was out of range on input.
c             return code i has precdence over i+1 for i = 9, 10, 11.
c iv(mlstgd) (i/o) saved value of iv(model).
c  iv(model) (i/o) on input, iv(model) should be an integer identifying
c             the current quadratic model of the objective function.
c             if a previous step yielded a better function reduction,
c             then iv(model) will be set to iv(mlstgd) on output.
c iv(nfcall) (in)  invocation count for the objective function.
c iv(nfgcal) (i/o) value of iv(nfcall) at step that gave the biggest
c             function reduction this iteration.  iv(nfgcal) remains
c             unchanged until a function reduction is obtained.
c iv(radinc) (i/o) the number of radius increases (or minus the number
c             of decreases) so far this iteration.
c iv(restor) (out) set to 0 unless x and v(f) have been restored, in
c             which case assess sets iv(restor) = 1.
c  iv(stage) (i/o) count of the number of models tried so far in the
c             current iteration.
c iv(stglim) (in)  maximum number of models to consider.
c iv(switch) (out) set to 0 unless a new model is being tried and it
c             gives a smaller function value than the previous model,
c             in which case assess sets iv(switch) = 1.
c iv(toobig) (in)  is nonzero if step was too big (e.g. if it caused
c             overflow).
c   iv(xirc) (i/o) value that iv(irc) would have in the absence of
c             convergence, false convergence, and oversized steps.
c
c  ***  v values referenced  ***
c
c v(afctol) (in)  absolute function convergence tolerance.  if the
c             absolute value of the current function value v(f) is less
c             than v(afctol), then assess returns with iv(irc) = 10.
c v(decfac) (in)  factor by which to decrease radius when iv(toobig) is
c             nonzero.
c v(dstnrm) (in)  the 2-norm of d*step.
c v(dstsav) (i/o) value of v(dstnrm) on saved step.
c   v(dst0) (in)  the 2-norm of d times the newton step (when defined,
c             i.e., for v(nreduc) .ge. 0).
c      v(f) (i/o) on both input and output, v(f) is the objective func-
c             tion value at x.  if x is restored to a previous value,
c             then v(f) is restored to the corresponding value.
c   v(fdif) (out) the function reduction v(f0) - v(f) (for the output
c             value of v(f) if an earlier step gave a bigger function
c             decrease, and for the input value of v(f) otherwise).
c v(flstgd) (i/o) saved value of v(f).
c     v(f0) (in)  objective function value at start of iteration.
c v(gtslst) (i/o) value of v(gtstep) on saved step.
c v(gtstep) (in)  inner product between step and gradient.
c v(incfac) (in)  minimum factor by which to increase radius.
c  v(lmax0) (in)  maximum reasonable step size (and initial step bound).
c             if the actual function decrease is no more than twice
c             what was predicted, if a return with iv(irc) = 7, 8, 9,
c             or 10 does not occur, if v(dstnrm) .gt. v(lmax0), and if
c             v(preduc) .le. v(rfctol) * abs(v(f0)), then assess re-
c             turns with iv(irc) = 11.  if so doing appears worthwhile,
c             then assess repeats this test with v(preduc) computed for
c             a step of length v(lmax0) (by a return with iv(irc) = 6).
c v(nreduc) (i/o)  function reduction predicted by quadratic model for
c             newton step.  if assess is called with iv(irc) = 6, i.e.,
c             if v(preduc) has been computed with radius = v(lmax0) for
c             use in the singular convervence test, then v(nreduc) is
c             set to -v(preduc) before the latter is restored.
c v(plstgd) (i/o) value of v(preduc) on saved step.
c v(preduc) (i/o) function reduction predicted by quadratic model for
c             current step.
c v(radfac) (out) factor to be used in determining the new radius,
c             which should be v(radfac)*dst, where  dst  is either the
c             output value of v(dstnrm) or the 2-norm of
c             diag(newd)*step  for the output value of step and the
c             updated version, newd, of the scale vector d.  for
c             iv(irc) = 3, v(radfac) = 1.0 is returned.
c v(rdfcmn) (in)  minimum value for v(radfac) in terms of the input
c             value of v(dstnrm) -- suggested value = 0.1.
c v(rdfcmx) (in)  maximum value for v(radfac) -- suggested value = 4.0.
c  v(reldx) (out) scaled relative change in x caused by step, computed
c             by function  reldst  as
c                 max (d(i)*abs(x(i)-x0(i)), 1 .le. i .le. p) /
c                    max (d(i)*(abs(x(i))+abs(x0(i))), 1 .le. i .le. p).
c             if an acceptable step is returned, then v(reldx) is com-
c             puted using the output (possibly restored) values of x
c             and step.  otherwise it is computed using the input
c             values.
c v(rfctol) (in)  relative function convergence tolerance.  if the
c             actual function reduction is at most twice what was pre-
c             dicted and  v(nreduc) .le. v(rfctol)*abs(v(f0)),  then
c             assess returns with iv(irc) = 8 or 9.  see also v(lmax0).
c v(stppar) (in)  marquardt parameter -- 0 means full newton step.
c v(tuner1) (in)  tuning constant used to decide if the function
c             reduction was much less than expected.  suggested
c             value = 0.1.
c v(tuner2) (in)  tuning constant used to decide if the function
c             reduction was large enough to accept step.  suggested
c             value = 10**-4.
c v(tuner3) (in)  tuning constant used to decide if the radius
c             should be increased.  suggested value = 0.75.
c  v(xctol) (in)  x-convergence criterion.  if step is a newton step
c             (v(stppar) = 0) having v(reldx) .le. v(xctol) and giving
c             at most twice the predicted function decrease, then
c             assess returns iv(irc) = 7 or 9.
c  v(xftol) (in)  false convergence tolerance.  if step gave no or only
c             a small function decrease and v(reldx) .le. v(xftol),
c             then assess returns with iv(irc) = 12.
c
c-------------------------------  notes  -------------------------------
c
c  ***  application and usage restrictions  ***
c
c        this routine is called as part of the nl2sol (nonlinear
c     least-squares) package.  it may be used in any unconstrained
c     minimization solver that uses dogleg, goldfeld-quandt-trotter,
c     or levenberg-marquardt steps.
c
c  ***  algorithm notes  ***
c
c        see (1) for further discussion of the assessing and model
c     switching strategies.  while nl2sol considers only two models,
c     assess is designed to handle any number of models.
c
c  ***  usage notes  ***
c
c        on the first call of an iteration, only the i/o variables
c     step, x, iv(irc), iv(model), v(f), v(dstnrm), v(gtstep), and
c     v(preduc) need have been initialized.  between calls, no i/o
c     values execpt step, x, iv(model), v(f) and the stopping toler-
c     ances should be changed.
c        after a return for convergence or false convergence, one can
c     change the stopping tolerances and call assess again, in which
c     case the stopping tests will be repeated.
c
c  ***  references  ***
c
c     (1) dennis, j.e., jr., gay, d.m., and welsch, r.e. (1981),
c        an adaptive nonlinear least-squares algorithm,
c        acm trans. math. software, vol. 7, no. 3.
c
c     (2) powell, m.j.d. (1970)  a fortran subroutine for solving
c        systems of nonlinear algebraic equations, in numerical
c        methods for nonlinear algebraic equations, edited by
c        p. rabinowitz, gordon and breach, london.
c
c  ***  history  ***
c
c        john dennis designed much of this routine, starting with
c     ideas in (2). roy welsch suggested the model switching strategy.
c        david gay and stephen peters cast this subroutine into a more
c     portable form (winter 1977), and david gay cast it into its
c     present form (fall 1978).
c
c  ***  general  ***
c
c     this subroutine was written in connection with research
c     supported by the national science foundation under grants
c     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and
c     mcs-7906671.
c
c------------------------  external quantities  ------------------------
c
c  ***  external functions and subroutines  ***
c
      external reldst, vcopy
      real reldst
c
c vcopy.... copies one vector to another.
c
c  ***  intrinsic functions  ***
c/+
      integer iabs
      real abs, amax1
c/
c  ***  no common blocks  ***
c
c--------------------------  local variables  --------------------------
c
      logical goodx
      integer i, nfc
      real emax, gts, half, one, reldx1, rfac1, two, xmax,
     1                 zero
c
c  ***  subscripts for iv and v  ***
c
      integer afctol, decfac, dstnrm, dstsav, dst0, f, fdif, flstgd, f0,
     1        gtslst, gtstep, incfac, irc, lmax0, mlstgd, model, nfcall,
     2        nfgcal, nreduc, plstgd, preduc, radfac, radinc, rdfcmn,
     3        rdfcmx, reldx, restor, rfctol, stage, stglim, stppar,
     4        switch, toobig, tuner1, tuner2, tuner3, xctol, xftol,
     5        xirc
c
c  ***  data initializations  ***
c
c/6
      data half/0.5e+0/, one/1.e+0/, two/2.e+0/, zero/0.e+0/
c/7
c     parameter (half=0.5d+0, one=1.d+0, two=2.d+0, zero=0.d+0)
c/
c
c/6
      data irc/3/, mlstgd/4/, model/5/, nfcall/6/,
     1     nfgcal/7/, radinc/8/, restor/9/, stage/10/,
     2     stglim/11/, switch/12/, toobig/2/, xirc/13/
c/7
c     parameter (irc=3, mlstgd=4, model=5, nfcall=6,
c    1     nfgcal=7, radinc=8, restor=9, stage=10,
c    2     stglim=11, switch=12, toobig=2, xirc=13)
c/
c/6
      data afctol/31/, decfac/22/, dstnrm/2/, dst0/3/,
     1     dstsav/18/, f/10/, fdif/11/, flstgd/12/, f0/13/,
     2     gtslst/14/, gtstep/4/, incfac/23/,
     3     lmax0/35/, nreduc/6/, plstgd/15/, preduc/7/,
     4     radfac/16/, rdfcmn/24/, rdfcmx/25/,
     5     reldx/17/, rfctol/32/, stppar/5/, tuner1/26/,
     6     tuner2/27/, tuner3/28/, xctol/33/, xftol/34/
c/7
c     parameter (afctol=31, decfac=22, dstnrm=2, dst0=3,
c    1     dstsav=18, f=10, fdif=11, flstgd=12, f0=13,
c    2     gtslst=14, gtstep=4, incfac=23,
c    3     lmax0=35, nreduc=6, plstgd=15, preduc=7,
c    4     radfac=16, rdfcmn=24, rdfcmx=25,
c    5     reldx=17, rfctol=32, stppar=5, tuner1=26,
c    6     tuner2=27, tuner3=28, xctol=33, xftol=34)
c/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
      nfc = iv(nfcall)
      iv(switch) = 0
      iv(restor) = 0
      rfac1 = one
      goodx = .true.
      i = iv(irc)
      if (i .ge. 1 .and. i .le. 12)
     1             go to (20,30,10,10,40,360,290,290,290,290,290,140), i
         iv(irc) = 13
         go to 999
c
c  ***  initialize for new iteration  ***
c
 10   iv(stage) = 1
      iv(radinc) = 0
      v(flstgd) = v(f0)
      if (iv(toobig) .eq. 0) go to 90
         iv(stage) = -1
         iv(xirc) = i
         go to 60
c
c  ***  step was recomputed with new model or smaller radius  ***
c  ***  first decide which  ***
c
 20   if (iv(model) .ne. iv(mlstgd)) go to 30
c        ***  old model retained, smaller radius tried  ***
c        ***  do not consider any more new models this iteration  ***
         iv(stage) = iv(stglim)
         iv(radinc) = -1
         go to 90
c
c  ***  a new model is being tried.  decide whether to keep it.  ***
c
 30   iv(stage) = iv(stage) + 1
c
c     ***  now we add the possibiltiy that step was recomputed with  ***
c     ***  the same model, perhaps because of an oversized step.     ***
c
 40   if (iv(stage) .gt. 0) go to 50
c
c        ***  step was recomputed because it was too big.  ***
c
         if (iv(toobig) .ne. 0) go to 60
c
c        ***  restore iv(stage) and pick up where we left off.  ***
c
         iv(stage) = -iv(stage)
         i = iv(xirc)
         go to (20, 30, 90, 90, 70), i
c
 50   if (iv(toobig) .eq. 0) go to 70
c
c  ***  handle oversize step  ***
c
      if (iv(radinc) .gt. 0) go to 80
         iv(stage) = -iv(stage)
         iv(xirc) = iv(irc)
c
 60      v(radfac) = v(decfac)
         iv(radinc) = iv(radinc) - 1
         iv(irc) = 5
         go to 999
c
 70   if (v(f) .lt. v(flstgd)) go to 90
c
c     *** the new step is a loser.  restore old model.  ***
c
      if (iv(model) .eq. iv(mlstgd)) go to 80
         iv(model) = iv(mlstgd)
         iv(switch) = 1
c
c     ***  restore step, etc. only if a previous step decreased v(f).
c
 80   if (v(flstgd) .ge. v(f0)) go to 90
         iv(restor) = 1
         v(f) = v(flstgd)
         v(preduc) = v(plstgd)
         v(gtstep) = v(gtslst)
         if (iv(switch) .eq. 0) rfac1 = v(dstnrm) / v(dstsav)
         v(dstnrm) = v(dstsav)
         nfc = iv(nfgcal)
         goodx = .false.
c
c
c  ***  compute relative change in x by current step  ***
c
 90   reldx1 = reldst(p, d, x, x0)
c
c  ***  restore x and step if necessary  ***
c
      if (goodx) go to 105
      do 100 i = 1, p
         step(i) = stlstg(i)
         x(i) = x0(i) + stlstg(i)
 100     continue
c
 105  v(fdif) = v(f0) - v(f)
      if (v(fdif) .gt. v(tuner2) * v(preduc)) go to 120
c
c        ***  no (or only a trivial) function decrease
c        ***  -- so try new model or smaller radius
c
         v(reldx) = reldx1
         if (v(f) .lt. v(f0)) go to 110
              iv(mlstgd) = iv(model)
              v(flstgd) = v(f)
              v(f) = v(f0)
              call vcopy(p, x, x0)
              iv(restor) = 1
              go to 115
 110     iv(nfgcal) = nfc
 115     iv(irc) = 1
         if (iv(stage) .lt. iv(stglim)) go to 130
              iv(irc) = 5
              iv(radinc) = iv(radinc) - 1
              go to 130
c
c  ***  nontrivial function decrease achieved  ***
c
 120  iv(nfgcal) = nfc
      rfac1 = one
      if (goodx) v(reldx) = reldx1
      v(dstsav) = v(dstnrm)
      if (v(fdif) .gt. v(preduc)*v(tuner1)) go to 200
c
c  ***  decrease was much less than predicted -- either change models
c  ***  or accept step with decreased radius.
c
      if (iv(stage) .ge. iv(stglim)) go to 125
c        ***  consider switching models  ***
         iv(irc) = 2
         go to 130
c
c     ***  accept step with decreased radius  ***
c
 125  iv(irc) = 4
c
c  ***  set v(radfac) to fletcher*s decrease factor  ***
c
 130  iv(xirc) = iv(irc)
      emax = v(gtstep) + v(fdif)
      v(radfac) = half * rfac1
      if (emax .lt. v(gtstep)) v(radfac) = rfac1 * amax1(v(rdfcmn),
     1                                           half * v(gtstep)/emax)
c
c  ***  do false convergence test  ***
c
 140  if (v(reldx) .le. v(xftol)) go to 160
         iv(irc) = iv(xirc)
         if (v(f) .lt. v(f0)) go to 230
              go to 300
c
 160  iv(irc) = 12
      go to 310
c
c  ***  handle good function decrease  ***
c
 200  if (v(fdif) .lt. (-v(tuner3) * v(gtstep))) go to 260
c
c     ***  increasing radius looks worthwhile.  see if we just
c     ***  recomputed step with a decreased radius or restored step
c     ***  after recomputing it with a larger radius.
c
      if (iv(radinc) .lt. 0) go to 260
      if (iv(restor) .eq. 1) go to 260
c
c        ***  we did not.  try a longer step unless this was a newton
c        ***  step.
c
         v(radfac) = v(rdfcmx)
         gts = v(gtstep)
         if (v(fdif) .lt. (half/v(radfac) - one) * gts)
     1            v(radfac) = amax1(v(incfac), half*gts/(gts + v(fdif)))
         iv(irc) = 4
         if (v(stppar) .eq. zero) go to 300
c             ***  step was not a newton step.  recompute it with
c             ***  a larger radius.
              iv(irc) = 5
              iv(radinc) = iv(radinc) + 1
c
c  ***  save values corresponding to good step  ***
c
 230  v(flstgd) = v(f)
      iv(mlstgd) = iv(model)
      call vcopy(p, stlstg, step)
      v(dstsav) = v(dstnrm)
      iv(nfgcal) = nfc
      v(plstgd) = v(preduc)
      v(gtslst) = v(gtstep)
      go to 300
c
c  ***  accept step with radius unchanged  ***
c
 260  v(radfac) = one
      iv(irc) = 3
      go to 300
c
c  ***  come here for a restart after convergence  ***
c
 290  iv(irc) = iv(xirc)
      if (v(dstsav) .ge. zero) go to 310
         iv(irc) = 12
         go to 310
c
c  ***  perform convergence tests  ***
c
 300  iv(xirc) = iv(irc)
 310  if (abs(v(f)) .lt. v(afctol)) iv(irc) = 10
      if (half * v(fdif) .gt. v(preduc)) go to 999
      emax = v(rfctol) * abs(v(f0))
      if (v(dstnrm) .gt. v(lmax0) .and. v(preduc) .le. emax)
     1                       iv(irc) = 11
      if (v(dst0) .lt. zero) go to 320
      i = 0
      if ((v(nreduc) .gt. zero .and. v(nreduc) .le. emax) .or.
     1    (v(nreduc) .eq. zero. and. v(preduc) .eq. zero))  i = 2
      if (v(stppar) .eq. zero .and. v(reldx) .le. v(xctol)
     1                        .and. goodx)                  i = i + 1
      if (i .gt. 0) iv(irc) = i + 6
c
c  ***  consider recomputing step of length v(lmax0) for singular
c  ***  convergence test.
c
 320  if (iabs(iv(irc)-3) .gt. 2 .and. iv(irc) .ne. 12) go to 999
      if (v(dstnrm) .gt. v(lmax0)) go to 330
         if (v(preduc) .ge. emax) go to 999
              if (v(dst0) .le. zero) go to 340
                   if (half * v(dst0) .le. v(lmax0)) go to 999
                        go to 340
 330  if (half * v(dstnrm) .le. v(lmax0)) go to 999
      xmax = v(lmax0) / v(dstnrm)
      if (xmax * (two - xmax) * v(preduc) .ge. emax) go to 999
 340  if (v(nreduc) .lt. zero) go to 370
c
c  ***  recompute v(preduc) for use in singular convergence test  ***
c
      v(gtslst) = v(gtstep)
      v(dstsav) = v(dstnrm)
      if (iv(irc) .eq. 12) v(dstsav) = -v(dstsav)
      v(plstgd) = v(preduc)
      iv(irc) = 6
      call vcopy(p, stlstg, step)
      go to 999
c
c  ***  perform singular convergence test with recomputed v(preduc)  ***
c
 360  v(gtstep) = v(gtslst)
      v(dstnrm) = abs(v(dstsav))
      call vcopy(p, step, stlstg)
      iv(irc) = iv(xirc)
      if (v(dstsav) .le. zero) iv(irc) = 12
      v(nreduc) = -v(preduc)
      v(preduc) = v(plstgd)
 370  if (-v(nreduc) .le. v(rfctol) * abs(v(f0))) iv(irc) = 11
c
 999  return
c
c  ***  last card of assess follows  ***
      end
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
      real d(p), j(nn,p), r(n), v(1), x(p)
c     dimension iv(*), v(*)
c
c  ***  local variables  ***
c
      logical havej
      integer cov, gp, gsave1, g1, hc, hmi, hpi, hpm, i, ipivi, ipivk,
     1        ip1, irc, k, kind, kl, l, m, mm1, mm1o2, pp1o2, qtr1,
     2        rd1, stpi, stpm, stp0, wl, w0, w1
      real del, half, negpt5, one, t, two, wk, zero
c
c  ***  intrinsic functions  ***
c/+
      integer iabs, max0
      real abs, amax1, float, sqrt
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
      data half/0.5e+0/, negpt5/-0.5e+0/, one/1.e+0/, two/2.e+0/,
     1     zero/0.e+0/
c/7
c     parameter (half=0.5d+0, negpt5=-0.5d+0, one=1.d+0, two=2.d+0,
c    1     zero=0.d+0)
c/
c
c/6
      data covmat/26/, covreq/15/, delta/50/, delta0/44/,
     1     dltfdc/40/, f/10/, fx/46/, g/28/, h/44/, ierr/32/,
     2     ipivot/61/, ipiv0/60/, kagqt/35/, kalm/36/,
     3     lmat/58/, mode/38/, nfgcal/7/, qtr/49/,
     4     rd/51/, rsave/52/, savei/54/, switch/12/,
     5     toobig/2/, w/59/, xmsave/49/
c/7
c     parameter (covmat=26, covreq=15, delta=50, delta0=44,
c    1     dltfdc=40, f=10, fx=46, g=28, h=44, ierr=32,
c    2     ipivot=61, ipiv0=60, kagqt=35, kalm=36,
c    3     lmat=58, mode=38, nfgcal=7, qtr=49,
c    4     rd=51, rsave=52, savei=54, switch=12,
c    5     toobig=2, w=59, xmsave=49)
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
      del = v(delta0) * amax1(one/d(m), abs(x(m)))
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
      del = v(dltfdc) * amax1(one/d(m), abs(x(m)))
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
      subroutine dfault(iv, v)
c
c  ***  supply nl2sol (version 2.2) default values to iv and v  ***
c
      integer iv(25)
      real v(45)
c/+
      real amax1
c/
      external imdcon, rmdcon
      integer imdcon
      real rmdcon
c
      real machep, mepcrt, one, sqteps, three
c
c  ***  subscripts for iv and v  ***
c
      integer afctol, cosmin, covprt, covreq, decfac, delta0, dfac,
     1        dinit, dltfdc, dltfdj, dtype, d0init, epslon, fuzz,
     2        incfac, inits, jtinit, lmax0, mxfcal, mxiter, outlev,
     3        parprt, phmnfc, phmxfc, prunit, rdfcmn, rdfcmx,
     4        rfctol, rlimit, solprt, statpr, tuner1, tuner2, tuner3,
     5        tuner4, tuner5, xctol, xftol, x0prt
c
c/6
      data one/1.e+0/, three/3.e+0/
c/7
c     parameter (one=1.d+0, three=3.d+0)
c/
c
c  ***  iv subscript values  ***
c
c/6
      data covprt/14/, covreq/15/, dtype/16/, inits/25/,
     1     mxfcal/17/, mxiter/18/, outlev/19/,
     2     parprt/20/, prunit/21/, solprt/22/,
     3     statpr/23/, x0prt/24/
c/7
c     parameter (covprt=14, covreq=15, dtype=16, inits=25,
c    1     mxfcal=17, mxiter=18, outlev=19,
c    2     parprt=20, prunit=21, solprt=22,
c    3     statpr=23, x0prt=24)
c/
c
c  ***  v subscript values  ***
c
c/6
      data afctol/31/, cosmin/43/, decfac/22/,
     1     delta0/44/, dfac/41/, dinit/38/, dltfdc/40/,
     2     dltfdj/36/, d0init/37/, epslon/19/, fuzz/45/,
     3     incfac/23/, jtinit/39/, lmax0/35/, phmnfc/20/,
     4     phmxfc/21/, rdfcmn/24/, rdfcmx/25/,
     5     rfctol/32/, rlimit/42/, tuner1/26/,
     6     tuner2/27/, tuner3/28/, tuner4/29/,
     7     tuner5/30/, xctol/33/, xftol/34/
c/7
c     parameter (afctol=31, cosmin=43, decfac=22,
c    1     delta0=44, dfac=41, dinit=38, dltfdc=40,
c    2     dltfdj=36, d0init=37, epslon=19, fuzz=45,
c    3     incfac=23, jtinit=39, lmax0=35, phmnfc=20,
c    4     phmxfc=21, rdfcmn=24, rdfcmx=25,
c    5     rfctol=32, rlimit=42, tuner1=26,
c    6     tuner2=27, tuner3=28, tuner4=29,
c    7     tuner5=30, xctol=33, xftol=34)
c/
c
c-----------------------------------------------------------------------
c
      iv(1) = 12
      iv(covprt) = 1
      iv(covreq) = 1
      iv(dtype) = 1
      iv(inits) = 0
      iv(mxfcal) = 200
      iv(mxiter) = 150
      iv(outlev) = 1
      iv(parprt) = 1
      iv(prunit) = imdcon(1)
      iv(solprt) = 1
      iv(statpr) = 1
      iv(x0prt) = 1
c
      machep = rmdcon(3)
      v(afctol) = 1.e-20
      if (machep .gt. 1.e-10) v(afctol) = machep**2
      v(cosmin) = amax1(1.e-6, 1.e+2 * machep)
      v(decfac) = 0.5e+0
      sqteps = rmdcon(4)
      v(delta0) = sqteps
      v(dfac) = 0.6e+0
      v(dinit) = 0.e+0
      mepcrt = machep ** (one/three)
      v(dltfdc) = mepcrt
      v(dltfdj) = sqteps
      v(d0init) = 1.e+0
      v(epslon) = 0.1e+0
      v(fuzz) = 1.5e+0
      v(incfac) = 2.e+0
      v(jtinit) = 1.e-6
      v(lmax0) = 100.e+0
      v(phmnfc) = -0.1e+0
      v(phmxfc) = 0.1e+0
      v(rdfcmn) = 0.1e+0
      v(rdfcmx) = 4.e+0
      v(rfctol) = amax1(1.e-10, mepcrt**2)
      v(rlimit) = rmdcon(5)
      v(tuner1) = 0.1e+0
      v(tuner2) = 1.e-4
      v(tuner3) = 0.75e+0
      v(tuner4) = 0.5e+0
      v(tuner5) = 0.75e+0
      v(xctol) = sqteps
      v(xftol) = 1.e+2 * machep
c
 999  return
c  ***  last card of dfault follows  ***
      end
      real function dotprd(p, x, y)
c
c  ***  return the inner product of the p-vectors x and y.  ***
c
      integer p
      real x(p), y(p)
c
      integer i
      real one, sqteta, t, zero
c/+
      real amax1, abs
c/
      external rmdcon
      real rmdcon
c
c  ***  rmdcon(2) returns a machine-dependent constant, sqteta, which
c  ***  is slightly larger than the smallest positive number that
c  ***  can be squared without underflowing.
c
c/6
      data one/1.e+0/, sqteta/0.e+0/, zero/0.e+0/
c/7
c     parameter (one=1.d+0, zero=0.d+0)
c     data sqteta/0.d+0/
c/
c
      dotprd = zero
      if (p .le. 0) go to 999
      if (sqteta .eq. zero) sqteta = rmdcon(2)
      do 20 i = 1, p
         t = amax1(abs(x(i)), abs(y(i)))
         if (t .gt. one) go to 10
         if (t .lt. sqteta) go to 20
         t = (x(i)/sqteta)*y(i)
         if (abs(t) .lt. sqteta) go to 20
 10      dotprd = dotprd + x(i)*y(i)
 20   continue
c
 999  return
c  ***  last card of dotprd follows  ***
      end
      subroutine dupdat(d, iv, j, n, nn, p, v)
c
c  ***  update scale vector d for nl2itr (nl2sol version 2.2)  ***
c
c  ***  parameter declarations  ***
c
      integer iv(1), n, nn, p
      real d(p), j(nn,p), v(1)
c     dimension iv(*), v(*)
c
c  ***  local variables  ***
c
      integer d0, i, jtoli, s1
      real sii, t, vdfac
c
c     ***  constants  ***
      real zero
c
c  ***  intrinsic functions  ***
c/+
      real amax1, sqrt
c/
c  ***  external function  ***
c
      external v2norm
      real v2norm
c
c  ***  subscripts for iv and v  ***
c
      integer dfac, dtype, jtol0, niter, s
c/6
      data dfac/41/, dtype/16/, jtol0/86/, niter/31/, s/53/
c/7
c     parameter (dfac=41, dtype=16, jtol0=86, niter=31, s=53)
c/
c
c/6
      data zero/0.e+0/
c/7
c     parameter (zero=0.d+0)
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
         if (sii .gt. zero) t = sqrt(t*t + sii)
         jtoli = jtol0 + i
         d0 = d0 + 1
         if (t .lt. v(jtoli)) t = amax1(v(d0), v(jtoli))
         d(i) = amax1(vdfac*d(i), t)
 30      continue
c
 999  return
c  ***  last card of dupdat follows  ***
      end
      subroutine gqtstp(d, dig, dihdi, ka, l, p, step, v, w)
c
c  *** compute goldfeld-quandt-trotter step by more-hebden technique ***
c  ***  (nl2sol version 2.2)  ***
c
c  ***  parameter declarations  ***
c
      integer ka, p
      real d(p), dig(p), dihdi(1), l(1), v(21), step(p),
     1                 w(1)
c     dimension dihdi(p*(p+1)/2), l(p*(p+1)/2), w(4*p+7)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  purpose  ***
c
c        given the (compactly stored) lower triangle of a scaled
c     hessian (approximation) and a nonzero scaled gradient vector,
c     this subroutine computes a goldfeld-quandt-trotter step of
c     approximate length v(radius) by the more-hebden technique.  in
c     other words, step is computed to (approximately) minimize
c     psi(step) = (g**t)*step + 0.5*(step**t)*h*step  such that the
c     2-norm of d*step is at most (approximately) v(radius), where
c     g  is the gradient,  h  is the hessian, and  d  is a diagonal
c     scale matrix whose diagonal is stored in the parameter d.
c     (gqtstp assumes  dig = d**-1 * g  and  dihdi = d**-1 * h * d**-1.)
c     if g = 0, however, step = 0 is returned (even at a saddle point).
c
c  ***  parameter description  ***
c
c     d (in)  = the scale vector, i.e. the diagonal of the scale
c              matrix  d  mentioned above under purpose.
c   dig (in)  = the scaled gradient vector, d**-1 * g.  if g = 0, then
c              step = 0  and  v(stppar) = 0  are returned.
c dihdi (in)  = lower triangle of the scaled hessian (approximation),
c              i.e., d**-1 * h * d**-1, stored compactly by rows., i.e.,
c              in the order (1,1), (2,1), (2,2), (3,1), (3,2), etc.
c    ka (i/o) = the number of hebden iterations (so far) taken to deter-
c              mine step.  ka .lt. 0 on input means this is the first
c              attempt to determine step (for the present dig and dihdi)
c              -- ka is initialized to 0 in this case.  output with
c              ka = 0  (or v(stppar) = 0)  means  step = -(h**-1)*g.
c     l (i/o) = workspace of length p*(p+1)/2 for cholesky factors.
c     p (in)  = number of parameters -- the hessian is a  p x p  matrix.
c  step (i/o) = the step computed.
c     v (i/o) contains various constants and variables described below.
c     w (i/o) = workspace of length 4*p + 6.
c
c  ***  entries in v  ***
c
c v(dgnorm) (i/o) = 2-norm of (d**-1)*g.
c v(dstnrm) (output) = 2-norm of d*step.
c v(dst0)   (i/o) = 2-norm of d*(h**-1)*g (for pos. def. h only), or
c             overestimate of smallest eigenvalue of (d**-1)*h*(d**-1).
c v(epslon) (in)  = max. rel. error allowed for psi(step).  for the
c             step returned, psi(step) will exceed its optimal value
c             by less than -v(epslon)*psi(step).  suggested value = 0.1.
c v(gtstep) (out) = inner product between g and step.
c v(nreduc) (out) = psi(-(h**-1)*g) = psi(newton step)  (for pos. def.
c             h only -- v(nreduc) is set to zero otherwise).
c v(phmnfc) (in)  = tol. (together with v(phmxfc)) for accepting step
c             (more*s sigma).  the error v(dstnrm) - v(radius) must lie
c             between v(phmnfc)*v(radius) and v(phmxfc)*v(radius).
c v(phmxfc) (in)  (see v(phmnfc).)
c             suggested values -- v(phmnfc) = -0.25, v(phmxfc) = 0.5.
c v(preduc) (out) = psi(step) = predicted obj. func. reduction for step.
c v(radius) (in)  = radius of current (scaled) trust region.
c v(rad0)   (i/o) = value of v(radius) from previous call.
c v(stppar) (i/o) is normally the marquardt parameter, i.e. the alpha
c             described below under algorithm notes.  if h + alpha*d**2
c             (see algorithm notes) is (nearly) singular, however,
c             then v(stppar) = -alpha.
c
c  ***  usage notes  ***
c
c     if it is desired to recompute step using a different value of
c     v(radius), then this routine may be restarted by calling it
c     with all parameters unchanged except v(radius).  (this explains
c     why step and w are listed as i/o).  on an intiial call (one with
c     ka .lt. 0), step and w need not be initialized and only compo-
c     nents v(epslon), v(stppar), v(phmnfc), v(phmxfc), v(radius), and
c     v(rad0) of v must be initialized.  to compute step from a saddle
c     point (where the true gradient vanishes and h has a negative
c     eigenvalue), a nonzero g with small components should be passed.
c
c  ***  application and usage restrictions  ***
c
c     this routine is called as part of the nl2sol (nonlinear least-
c     squares) package (ref. 1), but it could be used in solving any
c     unconstrained minimization problem.
c
c  ***  algorithm notes  ***
c
c        the desired g-q-t step (ref. 2, 3, 4) satisfies
c     (h + alpha*d**2)*step = -g  for some nonnegative alpha such that
c     h + alpha*d**2 is positive semidefinite.  alpha and step are
c     computed by a scheme analogous to the one described in ref. 5.
c     estimates of the smallest and largest eigenvalues of the hessian
c     are obtained from the gerschgorin circle theorem enhanced by a
c     simple form of the scaling described in ref. 6.  cases in which
c     h + alpha*d**2 is nearly (or exactly) singular are handled by
c     the technique discussed in ref. 2.  in these cases, a step of
c     (exact) length v(radius) is returned for which psi(step) exceeds
c     its optimal value by less than -v(epslon)*psi(step).
c
c  ***  functions and subroutines called  ***
c
c dotprd - returns inner product of two vectors.
c litvmu - applies inverse-transpose of compact lower triang. matrix.
c livmul - applies inverse of compact lower triang. matrix.
c lsqrt  - finds cholesky factor (of compactly stored lower triang.).
c lsvmin - returns approx. to min. sing. value of lower triang. matrix.
c rmdcon - returns machine-dependent constants.
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
c 3.  goldfeld, s.m., quandt, r.e., and trotter, h.f. (1966),
c             maximization by quadratic hill-climbing, econometrica 34,
c             pp. 541-551.
c 4.  hebden, m.d. (1973), an algorithm for minimization using exact
c             second derivatives, report t.p. 515, theoretical physics
c             div., a.e.r.e. harwell, oxon., england.
c 5.  more, j.j. (1978), the levenberg-marquardt algorithm, implemen-
c             tation and theory, pp.105-116 of springer lecture notes
c             in mathematics no. 630, edited by g.a. watson, springer-
c             verlag, berlin and new york.
c 6.  varga, r.s. (1965), minimal gerschgorin sets, pacific j. math. 15,
c             pp. 719-729.
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
      logical restrt
      integer dggdmx, diag, diag0, dstsav, emax, emin, i, im1, inc, irc,
     1        j, k, kalim, k1, lk0, phipin, q, q0, uk0, x, x0
      real alphak, aki, akk, delta, dst, epso6, lk,
     1                 oldphi, phi, phimax, phimin, psifac, rad,
     2                 root, si, sk, sw, t, twopsi, t1, uk, wi
c
c     ***  constants  ***
      real dgxfac, epsfac, four, half, kappa, negone, one,
     1                 p001, six, three, two, zero
c
c  ***  intrinsic functions  ***
c/+
      real abs, amax1, amin1, sqrt
c/
c  ***  external functions and subroutines  ***
c
      external dotprd, litvmu, livmul, lsqrt, lsvmin, rmdcon, v2norm
      real dotprd, lsvmin, rmdcon, v2norm
c
c  ***  subscripts for v  ***
c
      integer dgnorm, dstnrm, dst0, epslon, gtstep, stppar, nreduc,
     1        phmnfc, phmxfc, preduc, radius, rad0
c/6
      data dgnorm/1/, dstnrm/2/, dst0/3/, epslon/19/,
     1     gtstep/4/, nreduc/6/, phmnfc/20/,
     2     phmxfc/21/, preduc/7/, radius/8/,
     3     rad0/9/, stppar/5/
c/7
c     parameter (dgnorm=1, dstnrm=2, dst0=3, epslon=19,
c    1     gtstep=4, nreduc=6, phmnfc=20,
c    2     phmxfc=21, preduc=7, radius=8,
c    3     rad0=9, stppar=5)
c/
c
c/6
      data epsfac/50.0e+0/, four/4.0e+0/, half/0.5e+0/,
     1     kappa/2.0e+0/, negone/-1.0e+0/, one/1.0e+0/, p001/1.0e-3/,
     2     six/6.0e+0/, three/3.0e+0/, two/2.0e+0/, zero/0.0e+0/
c/7
c     parameter (epsfac=50.0d+0, four=4.0d+0, half=0.5d+0,
c    1     kappa=2.0d+0, negone=-1.0d+0, one=1.0d+0, p001=1.0d-3,
c    2     six=6.0d+0, three=3.0d+0, two=2.0d+0, zero=0.0d+0)
c     save dgxfac
c/
      data dgxfac/0.e+0/
c
c  ***  body  ***
c
c     ***  store largest abs. entry in (d**-1)*h*(d**-1) at w(dggdmx).
      dggdmx = p + 1
c     ***  store gerschgorin over- and underestimates of the largest
c     ***  and smallest eigenvalues of (d**-1)*h*(d**-1) at w(emax)
c     ***  and w(emin) respectively.
      emax = dggdmx + 1
      emin = emax + 1
c     ***  for use in recomputing step, the final values of lk, uk, dst,
c     ***  and the inverse derivative of more*s phi at 0 (for pos. def.
c     ***  h) are stored in w(lk0), w(uk0), w(dstsav), and w(phipin)
c     ***  respectively.
      lk0 = emin + 1
      phipin = lk0 + 1
      uk0 = phipin + 1
      dstsav = uk0 + 1
c     ***  store diag of (d**-1)*h*(d**-1) in w(diag),...,w(diag0+p).
      diag0 = dstsav
      diag = diag0 + 1
c     ***  store -d*step in w(q),...,w(q0+p).
      q0 = diag0 + p
      q = q0 + 1
      rad = v(radius)
c     ***  phitol = max. error allowed in dst = v(dstnrm) = 2-norm of
c     ***  d*step.
      phimax = v(phmxfc) * rad
      phimin = v(phmnfc) * rad
c     ***  epso6 and psifac are used in checking for the special case
c     ***  of (nearly) singular h + alpha*d**2 (see ref. 2).
      psifac = two * v(epslon) / (three * (four * (v(phmnfc) + one) *
     1                       (kappa + one)  +  kappa  +  two) * rad**2)
c     ***  oldphi is used to detect limits of numerical accuracy.  if
c     ***  we recompute step and it does not change, then we accept it.
      oldphi = zero
      epso6 = v(epslon)/six
      irc = 0
      restrt = .false.
      kalim = ka + 50
c
c  ***  start or restart, depending on ka  ***
c
      if (ka .ge. 0) go to 310
c
c  ***  fresh start  ***
c
      k = 0
      uk = negone
      ka = 0
      kalim = 50
c
c     ***  store diag(dihdi) in w(diag0+1),...,w(diag0+p)  ***
c
      j = 0
      do 20 i = 1, p
         j = j + i
         k1 = diag0 + i
         w(k1) = dihdi(j)
 20      continue
c
c     ***  determine w(dggdmx), the largest element of dihdi  ***
c
      t1 = zero
      j = p * (p + 1) / 2
      do 30 i = 1, j
         t = abs(dihdi(i))
         if (t1 .lt. t) t1 = t
 30      continue
      w(dggdmx) = t1
c
c  ***  try alpha = 0  ***
c
 40   call lsqrt(1, p, l, dihdi, irc)
      if (irc .eq. 0) go to 60
c        ***  indef. h -- underestimate smallest eigenvalue, use this
c        ***  estimate to initialize lower bound lk on alpha.
         j = irc*(irc+1)/2
         t = l(j)
         l(j) = one
         do 50 i = 1, irc
 50           w(i) = zero
         w(irc) = one
         call litvmu(irc, w, l, w)
         t1 = v2norm(irc, w)
         lk = -t / t1 / t1
         v(dst0) = -lk
         if (restrt) go to 210
         v(nreduc) = zero
         go to 70
c
c     ***  positive definite h -- compute unmodified newton step.  ***
 60   lk = zero
      call livmul(p, w(q), l, dig)
      v(nreduc) = half * dotprd(p, w(q), w(q))
      call litvmu(p, w(q), l, w(q))
      dst = v2norm(p, w(q))
      v(dst0) = dst
      phi = dst - rad
      if (phi .le. phimax) go to 280
      if (restrt) go to 210
c
c  ***  prepare to compute gerschgorin estimates of largest (and
c  ***  smallest) eigenvalues.  ***
c
 70   v(dgnorm) = v2norm(p, dig)
      if (v(dgnorm) .eq. zero) go to 450
      k = 0
      do 100 i = 1, p
         wi = zero
         if (i .eq. 1) go to 90
         im1 = i - 1
         do 80 j = 1, im1
              k = k + 1
              t = abs(dihdi(k))
              wi = wi + t
              w(j) = w(j) + t
 80           continue
 90      w(i) = wi
         k = k + 1
 100     continue
c
c  ***  (under-)estimate smallest eigenvalue of (d**-1)*h*(d**-1)  ***
c
      k = 1
      t1 = w(diag) - w(1)
      if (p .le. 1) go to 120
      do 110 i = 2, p
         j = diag0 + i
         t = w(j) - w(i)
         if (t .ge. t1) go to 110
              t1 = t
              k = i
 110     continue
c
 120  sk = w(k)
      j = diag0 + k
      akk = w(j)
      k1 = k*(k-1)/2 + 1
      inc = 1
      t = zero
      do 150 i = 1, p
         if (i .eq. k) go to 130
         aki = abs(dihdi(k1))
         si = w(i)
         j = diag0 + i
         t1 = half * (akk - w(j) + si - aki)
         t1 = t1 + sqrt(t1*t1 + sk*aki)
         if (t .lt. t1) t = t1
         if (i .lt. k) go to 140
 130     inc = i
 140     k1 = k1 + inc
 150     continue
c
      w(emin) = akk - t
      uk = v(dgnorm)/rad - w(emin)
c
c  ***  compute gerschgorin (over-)estimate of largest eigenvalue  ***
c
      k = 1
      t1 = w(diag) + w(1)
      if (p .le. 1) go to 170
      do 160 i = 2, p
         j = diag0 + i
         t = w(j) + w(i)
         if (t .le. t1) go to 160
              t1 = t
              k = i
 160     continue
c
 170  sk = w(k)
      j = diag0 + k
      akk = w(j)
      k1 = k*(k-1)/2 + 1
      inc = 1
      t = zero
      do 200 i = 1, p
         if (i .eq. k) go to 180
         aki = abs(dihdi(k1))
         si = w(i)
         j = diag0 + i
         t1 = half * (w(j) + si - aki - akk)
         t1 = t1 + sqrt(t1*t1 + sk*aki)
         if (t .lt. t1) t = t1
         if (i .lt. k) go to 190
 180     inc = i
 190     k1 = k1 + inc
 200     continue
c
      w(emax) = akk + t
      lk = amax1(lk, v(dgnorm)/rad - w(emax))
c
c     ***  alphak = current value of alpha (see alg. notes above).  we
c     ***  use more*s scheme for initializing it.
      alphak = abs(v(stppar)) * v(rad0)/rad
c
      if (irc .ne. 0) go to 210
c
c  ***  compute l0 for positive definite h  ***
c
      call livmul(p, w, l, w(q))
      t = v2norm(p, w)
      w(phipin) = dst / t / t
      lk = amax1(lk, phi*w(phipin))
c
c  ***  safeguard alphak and add alphak*i to (d**-1)*h*(d**-1)  ***
c
 210  ka = ka + 1
      if (-v(dst0) .ge. alphak .or. alphak .lt. lk .or. alphak .ge. uk)
     1                      alphak = uk * amax1(p001, sqrt(lk/uk))
      k = 0
      do 220 i = 1, p
         k = k + i
         j = diag0 + i
         dihdi(k) = w(j) + alphak
 220     continue
c
c  ***  try computing cholesky decomposition  ***
c
      call lsqrt(1, p, l, dihdi, irc)
      if (irc .eq. 0) go to 250
c
c  ***  (d**-1)*h*(d**-1) + alphak*i  is indefinite -- overestimate
c  ***  smallest eigenvalue for use in updating lk  ***
c
      j = (irc*(irc+1))/2
      t = l(j)
      l(j) = one
      do 230 i = 1, irc
 230     w(i) = zero
      w(irc) = one
      call litvmu(irc, w, l, w)
      t1 = v2norm(irc, w)
      lk = alphak - t/t1/t1
      v(dst0) = -lk
      go to 210
c
c  ***  alphak makes (d**-1)*h*(d**-1) positive definite.
c  ***  compute q = -d*step, check for convergence.  ***
c
 250  call livmul(p, w(q), l, dig)
      call litvmu(p, w(q), l, w(q))
      dst = v2norm(p, w(q))
      phi = dst - rad
      if (phi .le. phimax .and. phi .ge. phimin) go to 290
      if (phi .eq. oldphi) go to 290
      oldphi = phi
      if (phi .gt. zero) go to 260
c        ***  check for the special case of  h + alpha*d**2  (nearly)
c        ***  singular.  delta is .ge. the smallest eigenvalue of
c        ***  (d**-1)*h*(d**-1) + alphak*i.
         if (v(dst0) .gt. zero) go to 260
         delta = alphak + v(dst0)
         twopsi = alphak*dst*dst + dotprd(p, dig, w(q))
         if (delta .lt. psifac*twopsi) go to 270
c
c  ***  unacceptable alphak -- update lk, uk, alphak  ***
c
 260  if (ka .ge. kalim) go to 290
      call livmul(p, w, l, w(q))
      t1 = v2norm(p, w)
c     ***  the following dmin1 is necessary because of restarts  ***
      if (phi .lt. zero) uk = amin1(uk, alphak)
      alphak = alphak  +  (phi/t1) * (dst/t1) * (dst/rad)
      lk = amax1(lk, alphak)
      go to 210
c
c  ***  decide how to handle (nearly) singular h + alpha*d**2  ***
c
c     ***  if not yet available, obtain machine dependent value dgxfac.
 270  if (dgxfac .eq. zero) dgxfac = epsfac * rmdcon(3)
c
c     ***  now decide.  ***
      if (delta .gt. dgxfac*w(dggdmx)) go to 350
c        ***  delta is so small we cannot handle the special case in
c        ***  the available arithmetic.  accept step as it is.
         go to 290
c
c  ***  acceptable step on first try  ***
c
 280  alphak = zero
c
c  ***  successful step in general.  compute step = -(d**-1)*q  ***
c
 290  do 300 i = 1, p
         j = q0 + i
         step(i) = -w(j)/d(i)
 300     continue
      v(gtstep) = -dotprd(p, dig, w(q))
      v(preduc) = half * (abs(alphak)*dst*dst - v(gtstep))
      go to 430
c
c
c  ***  restart with new radius  ***
c
 310  if (v(dst0) .le. zero .or. v(dst0) - rad .gt. phimax) go to 330
c
c     ***  prepare to return newton step  ***
c
         restrt = .true.
         ka = ka + 1
         k = 0
         do 320 i = 1, p
              k = k + i
              j = diag0 + i
              dihdi(k) = w(j)
 320          continue
         uk = negone
         go to 40
c
 330  if (ka .eq. 0) go to 60
c
      dst = w(dstsav)
      alphak = abs(v(stppar))
      phi = dst - rad
      t = v(dgnorm)/rad
      if (rad .gt. v(rad0)) go to 340
c
c        ***  smaller radius  ***
         uk = t - w(emin)
         lk = zero
         if (alphak .gt. zero) lk = w(lk0)
         lk = amax1(lk, t - w(emax))
         if (v(dst0) .gt. zero) lk = amax1(lk, (v(dst0)-rad)*w(phipin))
         go to 260
c
c     ***  bigger radius  ***
 340  uk = t - w(emin)
      if (alphak .gt. zero) uk = amin1(uk, w(uk0))
      lk = amax1(zero, -v(dst0), t - w(emax))
      if (v(dst0) .gt. zero) lk = amax1(lk, (v(dst0)-rad)*w(phipin))
      go to 260
c
c  ***  handle (nearly) singular h + alpha*d**2  ***
c
c     ***  negate alphak to indicate special case  ***
 350  alphak = -alphak
c     ***  allocate storage for scratch vector x  ***
      x0 = q0 + p
      x = x0 + 1
c
c  ***  use inverse power method with start from lsvmin to obtain
c  ***  approximate eigenvector corresponding to smallest eigenvalue
c  ***  of (d**-1)*h*(d**-1).
c
      delta = kappa*delta
      t = lsvmin(p, l, w(x), w)
c
      k = 0
c     ***  normalize w  ***
 360  do 370 i = 1, p
 370     w(i) = t*w(i)
c     ***  complete current inv. power iter. -- replace w by (l**-t)*w.
      call litvmu(p, w, l, w)
      t1 = one/v2norm(p, w)
      t = t1*t
      if (t .le. delta) go to 390
      if (k .gt. 30) go to 290
      k = k + 1
c     ***  start next inv. power iter. by storing normalized w in x.
      do 380 i = 1, p
         j = x0 + i
         w(j) = t1*w(i)
 380     continue
c     ***  compute w = (l**-1)*x.
      call livmul(p, w, l, w(x))
      t = one/v2norm(p, w)
      go to 360
c
 390  do 400 i = 1, p
 400     w(i) = t1*w(i)
c
c  ***  now w is the desired approximate (unit) eigenvector and
c  ***  t*x = ((d**-1)*h*(d**-1) + alphak*i)*w.
c
      sw = dotprd(p, w(q), w)
      t1 = (rad + dst) * (rad - dst)
      root = sqrt(sw*sw + t1)
      if (sw .lt. zero) root = -root
      si = t1 / (sw + root)
c     ***  accept current step if adding si*w would lead to a
c     ***  further relative reduction in psi of less than v(epslon)/3.
      v(preduc) = half*twopsi
      t1 = zero
      t = si*(alphak*sw - half*si*(alphak + t*dotprd(p,w(x),w)))
      if (t .lt. epso6*twopsi) go to 410
         v(preduc) = v(preduc) + t
         dst = rad
         t1 = -si
 410  do 420 i = 1, p
         j = q0 + i
         w(j) = t1*w(i) - w(j)
         step(i) = w(j) / d(i)
 420     continue
      v(gtstep) = dotprd(p, dig, w(q))
c
c  ***  save values for use in a possible restart  ***
c
 430  v(dstnrm) = dst
      v(stppar) = alphak
      w(lk0) = lk
      w(uk0) = uk
      v(rad0) = rad
      w(dstsav) = dst
c
c     ***  restore diagonal of dihdi  ***
c
      j = 0
      do 440 i = 1, p
         j = j + i
         k = diag0 + i
         dihdi(j) = w(k)
 440     continue
      go to 999
c
c  ***  special case -- g = 0  ***
c
 450  v(stppar) = zero
      v(preduc) = zero
      v(dstnrm) = zero
      v(gtstep) = zero
      do 460 i = 1, p
 460     step(i) = zero
c
 999  return
c
c  ***  last card of gqtstp follows  ***
      end
      subroutine itsmry(d, iv, p, v, x)
c
c  ***  print nl2sol (version 2.2) iteration summary  ***
c
c  ***  parameter declarations  ***
c
      integer iv(1), p
      real d(p), v(1), x(p)
c     dimension iv(*), v(*)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  local variables  ***
c
      integer cov1, g1, i, ii, iv1, i1, j, m, nf, ng, ol, pu
c/6
      real model1(6), model2(6)
c/7
c     character*4 model1(6), model2(6)
c/
      real nreldf, oldf, preldf, reldf, zero
c
c  ***  intrinsic functions  ***
c/+
      integer iabs
c/
c  ***  no external functions or subroutines  ***
c
c  ***  subscripts for iv and v  ***
c
      integer covmat, covprt, covreq, dstnrm, f, fdif, f0, g,
     1        needhd, nfcall, nfcov, ngcov, ngcall, niter, nreduc,
     2        outlev, preduc, prntit, prunit, reldx, size, solprt,
     3        statpr, stppar, sused, x0prt
c
c  ***  iv subscript values  ***
c
c/6
      data covmat/26/, covprt/14/, g/28/, covreq/15/,
     1     needhd/39/, nfcall/6/, nfcov/40/, ngcov/41/,
     2     ngcall/30/, niter/31/, outlev/19/, prntit/48/,
     3     prunit/21/, solprt/22/, statpr/23/, sused/57/,
     4     x0prt/24/
c/7
c     parameter (covmat=26, covprt=14, g=28, covreq=15,
c    1     needhd=39, nfcall=6, nfcov=40, ngcov=41,
c    2     ngcall=30, niter=31, outlev=19, prntit=48,
c    3     prunit=21, solprt=22, statpr=23, sused=57,
c    4     x0prt=24)
c/
c
c  ***  v subscript values  ***
c
c/6
      data dstnrm/2/, f/10/, f0/13/, fdif/11/, nreduc/6/,
     1     preduc/7/, reldx/17/, size/47/, stppar/5/
c/7
c     parameter (dstnrm=2, f=10, f0=13, fdif=11, nreduc=6,
c    1     preduc=7, reldx=17, size=47, stppar=5)
c/
c
c/6
      data zero/0.e+0/
c/7
c     parameter (zero=0.d+0)
c/
c/6
      data model1(1)/4h    /, model1(2)/4h    /, model1(3)/4h    /,
     1     model1(4)/4h    /, model1(5)/4h  g /, model1(6)/4h  s /,
     2     model2(1)/4h g  /, model2(2)/4h s  /, model2(3)/4hg-s /,
     3     model2(4)/4hs-g /, model2(5)/4h-s-g/, model2(6)/4h-g-s/
c/7
c     data model1/'    ','    ','    ','    ','  g ','  s '/,
c    1     model2/' g  ',' s  ','g-s ','s-g ','-s-g','-g-s'/
c/
c
c-----------------------------------------------------------------------
c
      pu = iv(prunit)
      if (pu .eq. 0) go to 999
      iv1 = iv(1)
      ol = iv(outlev)
      if (iv1 .lt. 2 .or. iv1 .gt. 15) go to 140
      if (ol .eq. 0) go to 20
      if (iv1 .ge. 12) go to 20
      if (iv1 .ge. 10 .and. iv(prntit) .eq. 0) go to 20
      if (iv1 .gt. 2) go to 10
         iv(prntit) = iv(prntit) + 1
         if (iv(prntit) .lt. iabs(ol)) go to 999
 10   nf = iv(nfcall) - iabs(iv(nfcov))
      iv(prntit) = 0
      reldf = zero
      preldf = zero
      oldf = v(f0)
      if (oldf .le. zero) go to 12
         reldf = v(fdif) / oldf
         preldf = v(preduc) / oldf
 12   if (ol .gt. 0) go to 15
c
c        ***  print short summary line  ***
c
         if (iv(needhd) .eq. 1) write(pu, 1010)
 1010 format(12h0   it    nf,6x,1hf,8x,5hreldf,6x,6hpreldf,5x,5hreldx)
         iv(needhd) = 0
         write(pu,1017) iv(niter), nf, v(f), reldf, preldf, v(reldx)
         go to 20
c
c     ***  print long summary line  ***
c
 15   if (iv(needhd) .eq. 1) write(pu,1015)
 1015 format(12h0   it    nf,6x,1hf,8x,5hreldf,6x,6hpreldf,5x,5hreldx,
     1       4x,15hmodel    stppar,6x,4hsize,6x,6hd*step,5x,7hnpreldf)
      iv(needhd) = 0
      m = iv(sused)
      nreldf = zero
      if (oldf .gt. zero) nreldf = v(nreduc) / oldf
      write(pu,1017) iv(niter), nf, v(f), reldf, preldf, v(reldx),
     1               model1(m), model2(m), v(stppar), v(size),
     2               v(dstnrm), nreldf
 1017 format(1x,i5,i6,4e11.3,a3,a4,4e11.3)
c
 20   go to (999,999,30,35,40,45,50,60,70,80,90,150,110,120,130), iv1
c
 30   write(pu,1030)
 1030 format(26h0***** x-convergence *****)
      go to 180
c
 35   write(pu,1035)
 1035 format(42h0***** relative function convergence *****)
      go to 180
c
 40   write(pu,1040)
 1040 format(49h0***** x- and relative function convergence *****)
      go to 180
c
 45   write(pu,1045)
 1045 format(42h0***** absolute function convergence *****)
      go to 180
c
 50   write(pu,1050)
 1050 format(33h0***** singular convergence *****)
      go to 180
c
 60   write(pu,1060)
 1060 format(30h0***** false convergence *****)
      go to 180
c
 70   write(pu,1070)
 1070 format(38h0***** function evaluation limit *****)
      go to 180
c
 80   write(pu,1080)
 1080 format(28h0***** iteration limit *****)
      go to 180
c
 90   write(pu,1090)
 1090 format(18h0***** stopx *****)
      go to 180
c
 110  write(pu,1100)
 1100 format(45h0***** initial sum of squares overflows *****)
c
      go to 150
c
 120  write(pu,1120)
 1120 format(37h0***** bad parameters to assess *****)
      go to 999
c
 130  write(pu,1130)
 1130 format(36h0***** j could not be computed *****)
      if (iv(niter) .gt. 0) go to 190
      go to 150
c
 140  write(pu,1140) iv1
 1140 format(14h0***** iv(1) =,i5,6h *****)
      go to 999
c
c  ***  initial call on itsmry  ***
c
 150  if (iv(x0prt) .ne. 0) write(pu,1150) (i, x(i), d(i), i = 1, p)
 1150 format(23h0    i     initial x(i),7x,4hd(i)//(1x,i5,e17.6,e14.3))
      if (iv1 .ge. 13) go to 999
      iv(needhd) = 0
      iv(prntit) = 0
      if (ol .eq. 0) go to 999
      if (ol .lt. 0) write(pu,1010)
      if (ol .gt. 0) write(pu,1015)
      write(pu,1160) v(f)
 1160 format(12h0    0     1,e11.3,11x,e11.3)
      go to 999
c
c  ***  print various information requested on solution  ***
c
 180  iv(needhd) = 1
      if (iv(statpr) .eq. 0) go to 190
         oldf = v(f0)
         preldf = zero
         nreldf = zero
         if (oldf .le. zero) go to 185
              preldf = v(preduc) / oldf
              nreldf = v(nreduc) / oldf
 185     nf = iv(nfcall) - iv(nfcov)
         ng = iv(ngcall) - iv(ngcov)
         write(pu,1180) v(f), v(reldx), nf, ng, preldf, nreldf
 1180 format(9h0function,e17.6,8h   reldx,e20.6/12h func. evals,
     1   i8,9x,11hgrad. evals,i8/7h preldf,e19.6,3x,7hnpreldf,e18.6)
c
         if (iv(nfcov) .gt. 0) write(pu,1185) iv(nfcov)
 1185    format(1h0,i4,34h extra func. evals for covariance.)
         if (iv(ngcov) .gt. 0) write(pu,1186) iv(ngcov)
 1186    format(1x,i4,34h extra grad. evals for covariance.)
c
 190  if (iv(solprt) .eq. 0) go to 210
         iv(needhd) = 1
         g1 = iv(g)
         write(pu,1190)
 1190 format(22h0    i      final x(i),8x,4hd(i),10x,4hg(i)/)
         do 200 i = 1, p
              write(pu,1200) i, x(i), d(i), v(g1)
              g1 = g1 + 1
 200          continue
 1200    format(1x,i5,e17.6,2e14.3)
c
 210  if (iv(covprt) .eq. 0) go to 999
      cov1 = iv(covmat)
      iv(needhd) = 1
      if (cov1) 220, 230, 240
 220  if (-1 .eq. cov1) write(pu,1220)
 1220 format(43h0++++++ indefinite covariance matrix ++++++)
      if (-2 .eq. cov1) write(pu,1225)
 1225 format(52h0++++++ oversize steps in computing covariance +++++)
      go to 999
c
 230  write(pu,1230)
 1230 format(45h0++++++ covariance matrix not computed ++++++)
      go to 999
c
 240  i = iabs(iv(covreq))
      if (i .le. 1) write(pu,1241)
 1241 format(48h0covariance = scale * h**-1 * (j**t * j) * h**-1/)
      if (i .eq. 2) write(pu,1242)
 1242 format(27h0covariance = scale * h**-1/)
      if (i .ge. 3) write(pu,1243)
 1243 format(36h0covariance = scale * (j**t * j)**-1/)
      ii = cov1 - 1
      if (ol .le. 0) go to 260
      do 250 i = 1, p
         i1 = ii + 1
         ii = ii + i
         write(pu,1250) i, (v(j), j = i1, ii)
 250     continue
 1250 format(4h row,i3,2x,9e12.4/(9x,9e12.4))
      go to 999
c
 260  do 270 i = 1, p
         i1 = ii + 1
         ii = ii + i
         write(pu,1270) i, (v(j), j = i1, ii)
 270     continue
 1270 format(4h row,i3,2x,5e12.4/(9x,5e12.4))
c
 999  return
c  ***  last card of itsmry follows  ***
      end
      subroutine linvrt(n, lin, l)
c
c  ***  compute  lin = l**-1,  both  n x n  lower triang. stored   ***
c  ***  compactly by rows.  lin and l may share the same storage.  ***
c
c  ***  parameters  ***
c
      integer n
      real l(1), lin(1)
c     dimension l(n*(n+1)/2), lin(n*(n+1)/2)
c
c  ***  local variables  ***
c
      integer i, ii, im1, jj, j0, j1, k, k0, np1
      real one, t, zero
c/6
      data one/1.e+0/, zero/0.e+0/
c/7
c     parameter (one=1.d+0, zero=0.d+0)
c/
c
c  ***  body  ***
c
      np1 = n + 1
      j0 = n*(np1)/2
      do 30 ii = 1, n
         i = np1 - ii
         lin(j0) = one/l(j0)
         if (i .le. 1) go to 999
         j1 = j0
         im1 = i - 1
         do 20 jj = 1, im1
              t = zero
              j0 = j1
              k0 = j1 - jj
              do 10 k = 1, jj
                   t = t - l(k0)*lin(j0)
                   j0 = j0 - 1
                   k0 = k0 + k - i
 10                continue
              lin(j0) = t/l(k0)
 20           continue
         j0 = j0 - 1
 30      continue
 999  return
c  ***  last card of linvrt follows  ***
      end
      subroutine litvmu(n, x, l, y)
c
c  ***  solve  (l**t)*x = y,  where  l  is an  n x n  lower triangular
c  ***  matrix stored compactly by rows.  x and y may occupy the same
c  ***  storage.  ***
c
      integer n
      real x(n), l(1), y(n)
      integer i, ii, ij, im1, i0, j, np1
      real xi, zero
c/6
      data zero/0.e+0/
c/7
c     parameter (zero=0.d+0)
c/
c
      do 10 i = 1, n
 10      x(i) = y(i)
      np1 = n + 1
      i0 = n*(n+1)/2
      do 30 ii = 1, n
         i = np1 - ii
         xi = x(i)/l(i0)
         x(i) = xi
         if (i .le. 1) go to 999
         i0 = i0 - i
         if (xi .eq. zero) go to 30
         im1 = i - 1
         do 20 j = 1, im1
              ij = i0 + j
              x(j) = x(j) - xi*l(ij)
 20           continue
 30      continue
 999  return
c  ***  last card of litvmu follows  ***
      end
      subroutine livmul(n, x, l, y)
c
c  ***  solve  l*x = y, where  l  is an  n x n  lower triangular
c  ***  matrix stored compactly by rows.  x and y may occupy the same
c  ***  storage.  ***
c
      integer n
      real x(n), l(1), y(n)
      external dotprd
      real dotprd
      integer i, j, k
      real t, zero
c/6
      data zero/0.e+0/
c/7
c     parameter (zero=0.d+0)
c/
c
      do 10 k = 1, n
         if (y(k) .ne. zero) go to 20
         x(k) = zero
 10      continue
      go to 999
 20   j = k*(k+1)/2
      x(k) = y(k) / l(j)
      if (k .ge. n) go to 999
      k = k + 1
      do 30 i = k, n
         t = dotprd(i-1, l(j+1), x)
         j = j + i
         x(i) = (y(i) - t)/l(j)
 30      continue
 999  return
c  ***  last card of livmul follows  ***
      end
      subroutine lmstep(d, g, ierr, ipivot, ka, p, qtr, r, step, v, w)
c
c  ***  compute levenberg-marquardt step using more-hebden technique  **
c  ***  nl2sol version 2.2.  ***
c
c  ***  parameter declarations  ***
c
      integer ierr, ka, p
      integer ipivot(p)
      real d(p), g(p), qtr(p), r(1), step(p), v(21), w(1)
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
c 2.  gay, d.m. (1981), computing optimal locally constrained
c             siam j. sci. statist. computing, vol. 2, no. 2,
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
      real a, adi, alphak, b, dfacsq, dst, dtol, d1, d2,
     1                 lk, oldphi, phi, phimax, phimin, psifac, rad,
     2                 si, sj, sqrtak, t, twopsi, uk, wl
c
c     ***  constants  ***
      real dfac, eight, half, negone, one, p001, three,
     1                 ttol, zero
c
c  ***  intrinsic functions  ***
c/+
      integer iabs
      real abs, amax1, amin1, sqrt
c/
c  ***  external functions and subroutines  ***
c
      external dotprd, litvmu, livmul, vcopy, v2norm
      real dotprd, v2norm
c
c  ***  subscripts for v  ***
c
      integer dgnorm, dstnrm, dst0, epslon, gtstep, nreduc, phmnfc,
     1        phmxfc, preduc, radius, rad0, stppar
c/6
      data dgnorm/1/, dstnrm/2/, dst0/3/, epslon/19/,
     1     gtstep/4/, nreduc/6/, phmnfc/20/,
     2     phmxfc/21/, preduc/7/, radius/8/,
     3     rad0/9/, stppar/5/
c/7
c     parameter (dgnorm=1, dstnrm=2, dst0=3, epslon=19,
c    1     gtstep=4, nreduc=6, phmnfc=20,
c    2     phmxfc=21, preduc=7, radius=8,
c    3     rad0=9, stppar=5)
c/
c
c/6
      data dfac/256.e+0/, eight/8.e+0/, half/0.5e+0/, negone/-1.e+0/,
     1     one/1.e+0/, p001/1.e-3/, three/3.e+0/, ttol/2.5e+0/,
     2     zero/0.e+0/
c/7
c     parameter (dfac=256.d+0, eight=8.d+0, half=0.5d+0, negone=-1.d+0,
c    1     one=1.d+0, p001=1.d-3, three=3.d+0, ttol=2.5d+0,
c    2     zero=0.d+0)
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
      alphak = abs(v(stppar)) * v(rad0)/rad
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
     1             alphak = uk * amax1(p001, sqrt(lk/uk))
      sqrtak = sqrt(alphak)
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
         if (adi .ge. abs(wl)) go to 150
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
 190          if (abs(si) .gt. abs(wl)) go to 220
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
 300  if (phi .lt. zero) uk = amin1(uk, alphak)
      go to 320
 310  if (phi .lt. zero) uk = alphak
 320  do 330 i = 1, p
         j1 = ipivot(i)
         k = res0 + i
         step(i) = d(j1) * (w(k)/dst)
 330     continue
      call livmul(p, step, w(rmat), step)
      do 340 i = 1, p
 340     step(i) = step(i) / sqrt(w(i))
      t = one / v2norm(p, step)
      alphak = alphak + t*phi*t/rad
      lk = amax1(lk, alphak)
      go to 110
c
c  ***  restart  ***
c
 370  lk = w(lk0)
      uk = w(uk0)
      if (v(dst0) .gt. zero .and. v(dst0) - rad .le. phimax) go to 20
      alphak = abs(v(stppar))
      dst = w(dstsav)
      phi = dst - rad
      t = v(dgnorm)/rad
      if (rad .gt. v(rad0)) go to 380
c
c        ***  smaller radius  ***
         uk = t
         if (alphak .le. zero) lk = zero
         if (v(dst0) .gt. zero) lk = amax1(lk, (v(dst0)-rad)*w(phipin))
         go to 300
c
c     ***  bigger radius  ***
 380  if (alphak .le. zero .or. uk .gt. t) uk = t
      lk = zero
      if (v(dst0) .gt. zero) lk = amax1(lk, (v(dst0)-rad)*w(phipin))
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
      subroutine lsqrt(n1, n, l, a, irc)
c
c  ***  compute rows n1 through n of the cholesky factor  l  of
c  ***  a = l*(l**t),  where  l  and the lower triangle of  a  are both
c  ***  stored compactly by rows (and may occupy the same storage).
c  ***  irc = 0 means all went well.  irc = j means the leading
c  ***  principal  j x j  submatrix of  a  is not positive definite --
c  ***  and  l(j*(j+1)/2)  contains the (nonpos.) reduced j-th diagonal.
c
c  ***  parameters  ***
c
      integer n1, n, irc
      real l(1), a(1)
c     dimension l(n*(n+1)/2), a(n*(n+1)/2)
c
c  ***  local variables  ***
c
      integer i, ij, ik, im1, i0, j, jk, jm1, j0, k
      real t, td, zero
c
c  ***  intrinsic functions  ***
c/+
      real sqrt
c/
c/6
      data zero/0.e+0/
c/7
c     parameter (zero=0.d+0)
c/
c
c  ***  body  ***
c
      i0 = n1 * (n1 - 1) / 2
      do 50 i = n1, n
         td = zero
         if (i .eq. 1) go to 40
         j0 = 0
         im1 = i - 1
         do 30 j = 1, im1
              t = zero
              if (j .eq. 1) go to 20
              jm1 = j - 1
              do 10 k = 1, jm1
                   ik = i0 + k
                   jk = j0 + k
                   t = t + l(ik)*l(jk)
 10                continue
 20           ij = i0 + j
              j0 = j0 + j
              t = (a(ij) - t) / l(j0)
              l(ij) = t
              td = td + t*t
 30           continue
 40      i0 = i0 + i
         t = a(i0) - td
         if (t .le. zero) go to 60
         l(i0) = sqrt(t)
 50      continue
c
      irc = 0
      go to 999
c
 60   l(i0) = t
      irc = i
c
 999  return
c
c  ***  last card of lsqrt  ***
      end
      real function lsvmin(p, l, x, y)
c
c  ***  estimate smallest sing. value of packed lower triang. matrix l
c
c  ***  parameter declarations  ***
c
      integer p
      real l(1), x(p), y(p)
c     dimension l(p*(p+1)/2)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  purpose  ***
c
c     this function returns a good over-estimate of the smallest
c     singular value of the packed lower triangular matrix l.
c
c  ***  parameter description  ***
c
c  p (in)  = the order of l.  l is a  p x p  lower triangular matrix.
c  l (in)  = array holding the elements of  l  in row order, i.e.
c             l(1,1), l(2,1), l(2,2), l(3,1), l(3,2), l(3,3), etc.
c  x (out) if lsvmin returns a positive value, then x is a normalized
c             approximate left singular vector corresponding to the
c             smallest singular value.  this approximation may be very
c             crude.  if lsvmin returns zero, then some components of x
c             are zero and the rest retain their input values.
c  y (out) if lsvmin returns a positive value, then y = (l**-1)*x is an
c             unnormalized approximate right singular vector correspond-
c             ing to the smallest singular value.  this approximation
c             may be crude.  if lsvmin returns zero, then y retains its
c             input value.  the caller may pass the same vector for x
c             and y (nonstandard fortran usage), in which case y over-
c             writes x (for nonzero lsvmin returns).
c
c  ***  application and usage restrictions  ***
c
c     there are no usage restrictions.
c
c  ***  algorithm notes  ***
c
c     the algorithm is based on (1), with the additional provision that
c     lsvmin = 0 is returned if the smallest diagonal element of l
c     (in magnitude) is not more than the unit roundoff times the
c     largest.  the algorithm uses a random number generator proposed
c     in (4), which passes the spectral test with flying colors -- see
c     (2) and (3).
c
c  ***  subroutines and functions called  ***
c
c        v2norm - function, returns the 2-norm of a vector.
c
c  ***  references  ***
c
c     (1) cline, a., moler, c., stewart, g., and wilkinson, j.h.(1977),
c         an estimate for the condition number of a matrix, report
c         tm-310, applied math. div., argonne national laboratory.
c
c     (2) hoaglin, d.c. (1976), theoretical properties of congruential
c         random-number generators --  an empirical view,
c         memorandum ns-340, dept. of statistics, harvard univ.
c
c     (3) knuth, d.e. (1969), the art of computer programming, vol. 2
c         (seminumerical algorithms), addison-wesley, reading, mass.
c
c     (4) smith, c.s. (1971), multiplicative pseudo-random number
c         generators with prime modulus, j. assoc. comput. mach. 18,
c         pp. 586-593.
c
c  ***  history  ***
c
c     designed and coded by david m. gay (winter 1977/summer 1978).
c
c  ***  general  ***
c
c     this subroutine was written in connection with research
c     supported by the national science foundation under grants
c     mcs-7600324, dcr75-10143, 76-14311dss, and mcs76-11989.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  local variables  ***
c
      integer i, ii, ix, j, ji, jj, jjj, jm1, j0, pplus1
      real b, psj, sminus, splus, t, xminus, xplus
c
c  ***  constants  ***
c
      real half, one, r9973, zero
c
c  ***  intrinsic functions  ***
c/+
      integer mod
      real abs, float
c/
c  ***  external functions and subroutines  ***
c
      external v2norm
      real v2norm
c
c/6
      data half/0.5e+0/, one/1.e+0/, r9973/9973.e+0/, zero/0.e+0/
c/7
c     parameter (half=0.5d+0, one=1.d+0, r9973=9973.d+0, zero=0.d+0)
c     save ix
c/
      data ix/2/
c
c  ***  body  ***
c
c  ***  first check whether to return lsvmin = 0 and initialize x  ***
c
      ii = 0
      do 10 i = 1, p
         x(i) = zero
         ii = ii + i
         if (l(ii) .eq. zero) go to 300
 10      continue
      if (mod(ix, 9973) .eq. 0) ix = 2
      pplus1 = p + 1
c
c  ***  solve (l**t)*x = b, where the components of b have randomly
c  ***  chosen magnitudes in (.5,1) with signs chosen to make x large.
c
c     do j = p to 1 by -1...
      do 100 jjj = 1, p
         j = pplus1 - jjj
c       ***  determine x(j) in this iteration. note for i = 1,2,...,j
c       ***  that x(i) holds the current partial sum for row i.
         ix = mod(3432*ix, 9973)
         b = half*(one + float(ix)/r9973)
         xplus = (b - x(j))
         xminus = (-b - x(j))
         splus = abs(xplus)
         sminus = abs(xminus)
         jm1 = j - 1
         j0 = j*jm1/2
         jj = j0 + j
         xplus = xplus/l(jj)
         xminus = xminus/l(jj)
         if (jm1 .eq. 0) go to 30
         do 20 i = 1, jm1
              ji = j0 + i
              splus = splus + abs(x(i) + l(ji)*xplus)
              sminus = sminus + abs(x(i) + l(ji)*xminus)
 20           continue
 30      if (sminus .gt. splus) xplus = xminus
         x(j) = xplus
c       ***  update partial sums  ***
         if (jm1 .eq. 0) go to 100
         do 40 i = 1, jm1
              ji = j0 + i
              x(i) = x(i) + l(ji)*xplus
 40           continue
 100     continue
c
c  ***  normalize x  ***
c
      t = one/v2norm(p, x)
      do 110 i = 1, p
 110     x(i) = t*x(i)
c
c  ***  solve l*y = x and return svmin = 1/twonorm(y)  ***
c
      do 200 j = 1, p
         psj = zero
         jm1 = j - 1
         j0 = j*jm1/2
         if (jm1 .eq. 0) go to 130
         do 120 i = 1, jm1
              ji = j0 + i
              psj = psj + l(ji)*y(i)
 120          continue
 130     jj = j0 + j
         y(j) = (x(j) - psj)/l(jj)
 200     continue
c
      lsvmin = one/v2norm(p, y)
      go to 999
c
 300  lsvmin = zero
 999  return
c  ***  last card of lsvmin follows  ***
      end
      subroutine ltsqar(n, a, l)
c
c  ***  set a to lower triangle of (l**t) * l  ***
c
c  ***  l = n x n lower triang. matrix stored rowwise.  ***
c  ***  a is also stored rowwise and may share storage with l.  ***
c
      integer n
      real a(1), l(1)
c     dimension a(n*(n+1)/2), l(n*(n+1)/2)
c
      integer i, ii, iim1, i1, j, k, m
      real lii, lj
c
      ii = 0
      do 50 i = 1, n
         i1 = ii + 1
         ii = ii + i
         m = 1
         if (i .eq. 1) go to 30
         iim1 = ii - 1
         do 20 j = i1, iim1
              lj = l(j)
              do 10 k = i1, j
                   a(m) = a(m) + lj*l(k)
                   m = m + 1
 10                continue
 20           continue
 30      lii = l(ii)
         do 40 j = i1, ii
 40           a(j) = lii * l(j)
 50      continue
c
 999  return
c  ***  last card of ltsqar follows  ***
      end
      subroutine parchk(iv, n, nn, p, v)
c
c  ***  check nl2sol (version 2.2) parameters, print changed values  ***
c
      integer iv(1), n, nn, p
      real v(1)
c     dimension iv(*), v(*)
c
      external dfault, rmdcon, vcopy
      real rmdcon
c dfault -- supplies dfault parameter values.
c rmdcon -- returns machine-dependent constants.
c vcopy  -- copies one vector to another.
c
c  ***  local variables  ***
c
      integer i, iv1, jtolp, k, l, m, nvdflt, pu
c/6
      real cngd(3), dflt(3), vn(2,27), which(3)
c/7
c     character*4 cngd(3), dflt(3), vn(2,27), which(3)
c/
      real big, machep, tiny, vk, vm(27), vx(27), zero
c
c  ***  iv and v subscripts  ***
c
      integer dtype, dtype0, d0init, epslon, inits, jtinit, jtol0,
     1        jtol1, oldn, oldnn, oldp, parprt, parsv1, prunit
c
c/6
      data nvdflt/27/, zero/0.e+0/
c/7
c     parameter (nvdflt=27, zero=0.d+0)
c/
c
c/6
      data dtype/16/, dtype0/29/, d0init/37/, epslon/19/,
     1     inits/25/, jtinit/39/, jtol0/86/, jtol1/87/,
     2     oldn/45/, oldnn/46/, oldp/47/, parprt/20/,
     3     parsv1/51/, prunit/21/
c/7
c     parameter (dtype=16, dtype0=29, d0init=37, epslon=19,
c    1     inits=25, jtinit=39, jtol0=86, jtol1=87,
c    2     oldn=45, oldnn=46, oldp=47, parprt=20,
c    3     parsv1=51, prunit=21)
c     save big, tiny
c/
c
      data big/0.e+0/, tiny/1.e+0/
c/6
      data vn(1,1),vn(2,1)/4hepsl,4hon../
      data vn(1,2),vn(2,2)/4hphmn,4hfc../
      data vn(1,3),vn(2,3)/4hphmx,4hfc../
      data vn(1,4),vn(2,4)/4hdecf,4hac../
      data vn(1,5),vn(2,5)/4hincf,4hac../
      data vn(1,6),vn(2,6)/4hrdfc,4hmn../
      data vn(1,7),vn(2,7)/4hrdfc,4hmx../
      data vn(1,8),vn(2,8)/4htune,4hr1../
      data vn(1,9),vn(2,9)/4htune,4hr2../
      data vn(1,10),vn(2,10)/4htune,4hr3../
      data vn(1,11),vn(2,11)/4htune,4hr4../
      data vn(1,12),vn(2,12)/4htune,4hr5../
      data vn(1,13),vn(2,13)/4hafct,4hol../
      data vn(1,14),vn(2,14)/4hrfct,4hol../
      data vn(1,15),vn(2,15)/4hxcto,4hl.../
      data vn(1,16),vn(2,16)/4hxfto,4hl.../
      data vn(1,17),vn(2,17)/4hlmax,4h0.../
      data vn(1,18),vn(2,18)/4hdltf,4hdj../
      data vn(1,19),vn(2,19)/4hd0in,4hit../
      data vn(1,20),vn(2,20)/4hdini,4ht.../
      data vn(1,21),vn(2,21)/4hjtin,4hit../
      data vn(1,22),vn(2,22)/4hdltf,4hdc../
      data vn(1,23),vn(2,23)/4hdfac,4h..../
      data vn(1,24),vn(2,24)/4hrlim,4hit../
      data vn(1,25),vn(2,25)/4hcosm,4hin../
      data vn(1,26),vn(2,26)/4hdelt,4ha0../
      data vn(1,27),vn(2,27)/4hfuzz,4h..../
c/7
c     data vn(1,1),vn(2,1)/'epsl','on..'/
c     data vn(1,2),vn(2,2)/'phmn','fc..'/
c     data vn(1,3),vn(2,3)/'phmx','fc..'/
c     data vn(1,4),vn(2,4)/'decf','ac..'/
c     data vn(1,5),vn(2,5)/'incf','ac..'/
c     data vn(1,6),vn(2,6)/'rdfc','mn..'/
c     data vn(1,7),vn(2,7)/'rdfc','mx..'/
c     data vn(1,8),vn(2,8)/'tune','r1..'/
c     data vn(1,9),vn(2,9)/'tune','r2..'/
c     data vn(1,10),vn(2,10)/'tune','r3..'/
c     data vn(1,11),vn(2,11)/'tune','r4..'/
c     data vn(1,12),vn(2,12)/'tune','r5..'/
c     data vn(1,13),vn(2,13)/'afct','ol..'/
c     data vn(1,14),vn(2,14)/'rfct','ol..'/
c     data vn(1,15),vn(2,15)/'xcto','l...'/
c     data vn(1,16),vn(2,16)/'xfto','l...'/
c     data vn(1,17),vn(2,17)/'lmax','0...'/
c     data vn(1,18),vn(2,18)/'dltf','dj..'/
c     data vn(1,19),vn(2,19)/'d0in','it..'/
c     data vn(1,20),vn(2,20)/'dini','t...'/
c     data vn(1,21),vn(2,21)/'jtin','it..'/
c     data vn(1,22),vn(2,22)/'dltf','dc..'/
c     data vn(1,23),vn(2,23)/'dfac','....'/
c     data vn(1,24),vn(2,24)/'rlim','it..'/
c     data vn(1,25),vn(2,25)/'cosm','in..'/
c     data vn(1,26),vn(2,26)/'delt','a0..'/
c     data vn(1,27),vn(2,27)/'fuzz','....'/
c/
c
      data vm(1)/1.0e-3/, vm(2)/-0.99e+0/, vm(3)/1.0e-3/, vm(4)/1.0e-2/,
     1     vm(5)/1.2e+0/, vm(6)/1.e-2/, vm(7)/1.2e+0/, vm(8)/0.e+0/,
     2     vm(9)/0.e+0/, vm(10)/1.e-3/, vm(11)/-1.e+0/, vm(15)/0.e+0/,
     3     vm(16)/0.e+0/, vm(19)/0.e+0/, vm(20)/-10.e+0/, vm(21)/0.e+0/,
     4     vm(23)/0.e+0/, vm(24)/1.e+10/, vm(27)/1.01e+0/
      data vx(1)/0.9e+0/, vx(2)/-1.e-3/, vx(3)/1.e+1/, vx(4)/0.8e+0/,
     1     vx(5)/1.e+2/, vx(6)/0.8e+0/, vx(7)/1.e+2/, vx(8)/0.5e+0/,
     2     vx(9)/0.5e+0/, vx(10)/1.e+0/, vx(11)/1.e+0/, vx(14)/0.1e+0/,
     3     vx(15)/1.e+0/, vx(16)/1.e+0/, vx(18)/1.e+0/, vx(22)/1.e+0/,
     4     vx(23)/1.e+0/, vx(25)/1.e+0/, vx(26)/1.e+0/, vx(27)/1.e+2/
c
c/6
      data cngd(1),cngd(2),cngd(3)/4h---c,4hhang,4hed v/,
     1     dflt(1),dflt(2),dflt(3)/4hnond,4hefau,4hlt v/
c/7
c     data cngd(1),cngd(2),cngd(3)/'---c','hang','ed v'/,
c    1     dflt(1),dflt(2),dflt(3)/'nond','efau','lt v'/
c/
c
c.......................................................................
c
      if (iv(1) .eq. 0) call dfault(iv, v)
      pu = iv(prunit)
      iv1 = iv(1)
      if (iv1 .ne. 12) go to 30
         if (nn .ge. n .and. n .ge. p .and. p .ge. 1) go to 20
              iv(1) = 16
              if (pu .ne. 0) write(pu,10) nn, n, p
 10           format(30h0///// bad nn, n, or p... nn =,i5,5h, n =,i5,
     1               5h, p =,i5)
              go to 999
 20      k = iv(21)
         call dfault(iv(21), v(33))
         iv(21) = k
         iv(dtype0) = iv(dtype+20)
         iv(oldn) = n
         iv(oldnn) = nn
         iv(oldp) = p
         which(1) = dflt(1)
         which(2) = dflt(2)
         which(3) = dflt(3)
         go to 80
 30   if (n .eq. iv(oldn) .and. nn .eq. iv(oldnn) .and. p .eq. iv(oldp))
     1                       go to 50
         iv(1) = 17
         if (pu .ne. 0) write(pu,40) iv(oldnn), iv(oldn), iv(oldp), nn,
     1                               n, p
 40      format(30h0///// (nn,n,p) changed from (,i5,1h,,i5,1h,,i3,
     1          6h) to (,i5,1h,,i5,1h,,i3,2h).)
         go to 999
c
 50   if (iv1 .le. 11 .and. iv1 .ge. 1) go to 70
         iv(1) = 50
         if (pu .ne. 0) write(pu,60) iv1
 60      format(15h0/////  iv(1) =,i5,28h should be between 0 and 12.)
         go to 999
c
 70   which(1) = cngd(1)
      which(2) = cngd(2)
      which(3) = cngd(3)
c
 80   if (big .gt. tiny) go to 90
         tiny = rmdcon(1)
         machep = rmdcon(3)
         big = rmdcon(6)
         vm(12) = machep
         vx(12) = big
         vm(13) = tiny
         vx(13) = big
         vm(14) = machep
         vm(17) = tiny
         vx(17) = big
         vm(18) = machep
         vx(19) = big
         vx(20) = big
         vx(21) = big
         vm(22) = machep
         vx(24) = rmdcon(5)
         vm(25) = machep
         vm(26) = machep
 90   m = 0
      if (iv(inits) .ge. 0 .and. iv(inits) .le. 2) go to 110
         m = 18
         if (pu .ne. 0) write(pu,100) iv(inits)
 100     format(25h0/////  inits... iv(25) =,i4,20h should be between 0,
     1          7h and 2.)
 110  k = epslon
      do 140 i = 1, nvdflt
         vk = v(k)
         if (vk .ge. vm(i) .and. vk .le. vx(i)) go to 130
              m = k
              if (pu .ne. 0) write(pu,120) vn(1,i), vn(2,i), k, vk,
     1                                    vm(i), vx(i)
 120          format(8h0/////  ,2a4,5h.. v(,i2,3h) =,e11.3,7h should,
     1               11h be between,e11.3,4h and,d11.3)
 130     k = k + 1
 140     continue
c
      if (iv1 .eq. 12 .and. v(jtinit) .gt. zero) go to 170
c
c  ***  check jtol values  ***
c
      jtolp = jtol0 + p
      do 160 i = jtol1, jtolp
         if (v(i) .gt. zero) go to 160
         k = i - jtol0
         if (pu .ne. 0) write(pu,150) k, i, v(i)
 150     format(12h0///// jtol(,i3,6h) = v(,i3,3h) =,e11.3,
     1          20h should be positive.)
         m = i
 160     continue
c
 170  if (m .eq. 0) go to 180
         iv(1) = m
         go to 999
c
 180  if (pu .eq. 0 .or. iv(parprt) .eq. 0) go to 999
      if (iv1 .ne. 12 .or. iv(inits) .eq. 0) go to 200
         m = 1
         write(pu,190) iv(inits)
 190     format(22h0nondefault values..../20h inits..... iv(25) =,i3)
 200  if (iv(dtype) .eq. iv(dtype0)) go to  210
         if (m .eq. 0) write(pu,215) which
         m = 1
         write(pu,205) iv(dtype)
 205     format(20h dtype..... iv(16) =,i3)
 210  k = epslon
      l = parsv1
      do 240 i = 1, nvdflt
         if (v(k) .eq. v(l)) go to 230
              if (m .eq. 0) write(pu,215) which
 215          format(1h0,3a4,9halues..../)
              m = 1
              write(pu,220) vn(1,i), vn(2,i), k, v(k)
 220          format(1x,2a4,5h.. v(,i2,3h) =,e15.7)
 230     k = k + 1
         l = l + 1
 240     continue
      iv(dtype0) = iv(dtype)
      call vcopy(nvdflt, v(parsv1), v(epslon))
      if (iv1 .ne. 12) go to 999
         if (v(jtinit) .gt. zero) go to 260
              jtolp = jtol0 + p
              write(pu,250) (v(i), i = jtol1, jtolp)
 250          format(24h0(initial) jtol array.../(1x,6e12.3))
 260     if (v(d0init) .gt. zero) go to 999
              k = jtol1 + p
              l = k + p - 1
              write(pu,270) (v(i), i = k, l)
 270          format(22h0(initial) d0 array.../1x,6e12.3)
c
 999  return
c  ***  last card of parchk follows  ***
      end
      subroutine qapply(nn, n, p, j, r, ierr)
c     *****parameters.
      integer nn, n, p, ierr
      real j(nn,p), r(n)
c
c     ..................................................................
c     ..................................................................
c
c     *****purpose.
c     this subroutine applies to r the orthogonal transformations
c     stored in j by qrfact
c
c     *****parameter description.
c     on input.
c
c        nn is the row dimension of the matrix j as declared in
c             the calling program dimension statement
c
c        n is the number of rows of j and the size of the vector r
c
c        p is the number of columns of j and the size of sigma
c
c        j contains on and below its diagonal the column vectors
c             u which determine the householder transformations
c             ident - u*u.transpose
c
c        r is the right hand side vector to which the orthogonal
c             transformations will be applied
c
c        ierr if non-zero indicates that not all the transformations
c             were successfully determined and only the first
c             abs(ierr) - 1 transformations will be used
c
c     on output.
c
c        r has been overwritten by its transformed image
c
c     *****application and usage restrictions.
c     none
c
c     *****algorithm notes.
c     the vectors u which determine the householder transformations
c     are normalized so that their 2-norm squared is 2.  the use of
c     these transformations here is in the spirit of (1).
c
c     *****subroutines and functions called.
c
c     dotprd - function, returns the inner product of vectors
c
c     *****references.
c     (1) businger, p. a., and golub, g. h. (1965), linear least squares
c        solutions by householder transformations, numer. math. 7,
c        pp. 269-276.
c
c     *****history.
c     designed by david m. gay, coded by stephen c. peters (winter 1977)
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
c     *****local variables.
      integer i, k, l, nl1
      real t
c     *****intrinsic functions.
c/+
      integer iabs
c/
c     *****functions.
      external dotprd
      real dotprd
c
      k = p
      if (ierr .ne. 0) k = iabs(ierr) - 1
      if ( k .eq. 0) go to 999
c
      do 20 l = 1, k
         nl1 = n - l + 1
         t = -dotprd(nl1, j(l,l), r(l))
c
         do 10 i = l, n
 10           r(i) = r(i) + t*j(i,l)
 20   continue
 999  return
c     .... last card of qapply .........................................
      end
      subroutine qrfact(nm,m,n,qr,alpha,ipivot,ierr,nopivk,sum)
c
c  ***  compute the qr decomposition of the matrix stored in qr  ***
c
c     *****parameters.
      integer nm,m,n,ipivot(n),ierr,nopivk
      real              qr(nm,n),alpha(n),sum(n)
c     *****local variables.
      integer i,j,jbar,k,k1,minum,mk1
      real              alphak,beta,qrkk,qrkmax,sigma,temp,ufeta,rktol,
     1        rktol1,sumj
c     *****functions.
c/+
      integer min0
      real              abs,sqrt
c/
      external dotprd, rmdcon, vaxpy, vscopy, v2norm
      real dotprd, rmdcon, v2norm
c dotprd... returns inner product of two vectors.
c rmdcon... returns machine-dependent constants.
c vaxpy... computes scalar times one vector plus another.
c vscopy... sets all elements of a vector to a scalar.
c v2norm... returns the 2-norm of a vector.
c
c     *****constants.
      real one, p01, p99, zero
c/6
      data one/1.0e+0/, p01/0.01e+0/, p99/0.99e+0/, zero/0.0e+0/
c/7
c     parameter (one=1.0d+0, p01=0.01d+0, p99=0.99d+0, zero=0.0d+0)
c     save rktol, ufeta
c/
c
c
c     ..................................................................
c     ..................................................................
c
c
c     *****purpose.
c
c     this subroutine does a qr-decomposition on the m x n matrix qr,
c        with an optionally modified column pivoting, and returns the
c        upper triangular r-matrix, as well as the orthogonal vectors
c        used in the transformations.
c
c     *****parameter description.
c     on input.
c
c        nm must be set to the row dimension of the two dimensional
c             array parameters as declared in the calling program
c             dimension statement.
c
c        m must be set to the number of rows in the matrix.
c
c        n must be set to the number of columns in the matrix.
c
c        qr contains the real rectangular matrix to be decomposed.
c
c     nopivk is used to control pivotting.  columns 1 through
c        nopivk will remain fixed in position.
c
c        sum is used for temporary storage for the subroutine.
c
c     on output.
c
c        qr contains the non-diagonal elements of the r-matrix
c             in the strict upper triangle. the vectors u, which
c             define the householder transformations   i - u*u-transp,
c             are in the columns of the lower triangle. these vectors u
c             are scaled so that the square of their 2-norm is 2.0.
c
c        alpha contains the diagonal elements of the r-matrix.
c
c        ipivot reflects the column pivoting performed on the input
c             matrix to accomplish the decomposition. the j-th
c             element of ipivot gives the column of the original
c             matrix which was pivoted into column j during the
c             decomposition.
c
c        ierr is set to.
c             0 for normal return,
c             k if no non-zero pivot could be found for the k-th
c                  transformation, or
c             -k for an error exit on the k-th thansformation.
c             if an error exit was taken, the first (k - 1)
c             transformations are correct.
c
c
c     *****applications and usage restrictions.
c     this may be used when solving linear least-squares problems --
c     see subroutine qr1 of rosepack.  it is called for this purpose
c     by llsqst in the nl2sol (nonlinear least-squares) package.
c
c     *****algorithm notes.
c     this version of qrfact tries to eliminate the occurrence of
c     underflows during the accumulation of inner products.  rktol1
c     is chosen below so as to insure that discarded terms have no
c     effect on the computed two-norms.
c
c     adapted from the algol routine solve (1).
c
c     *****references.
c     (1)     businger,p. and golub,g.h., linear least squares
c     solutions by housholder transformations, in wilkinson,j.h.
c     and reinsch,c.(eds.), handbook for automatic computation,
c     volume ii. linear algebra, springer-verlag, 111-118 (1971).
c     prepublished in numer.math. 7, 269-276 (1965).
c
c     *****history.
c     this amounts to the subroutine qr1 of rosepack with rktol1 used
c     in place of rktol below, with v2norm used to initialize (and
c     sometimes update) the sum array, and with calls on dotprd and
c     vaxpy in place of some loops.
c
c     *****general.
c
c     development of this program supported in part by
c     national science foundation grant gj-1154x3 and
c     national science foundation grant dcr75-08802
c     to national bureau of economic research, inc.
c
c
c
c     ..................................................................
c     ..................................................................
c
c
c     ..........  ufeta is the smallest positive floating point number
c        s.t. ufeta and -ufeta can both be represented.
c
c     ..........  rktol is the square root of the relative precision
c        of floating point arithmetic (machep).
      data rktol/0.e+0/, ufeta/0.e+0/
c     *****body of program.
      if (ufeta .gt. zero) go to 10
         ufeta = rmdcon(1)
         rktol = rmdcon(4)
   10 ierr = 0
      rktol1 = p01 * rktol
c
      do 20 j=1,n
         sum(j) = v2norm(m, qr(1,j))
         ipivot(j) = j
   20 continue
c
      minum = min0(m,n)
c
      do 120 k=1,minum
         mk1 = m - k + 1
c        ..........k-th householder transformation..........
         sigma = zero
         jbar = 0
c        ..........find largest column sum..........
      if (k .le. nopivk) go to 50
         do 30 j=k,n
              if (sigma .ge. sum(j))  go to 30
              sigma = sum(j)
              jbar = j
   30    continue
c
         if (jbar .eq. 0)  go to 220
         if (jbar .eq. k)  go to 50
c        ..........column interchange..........
         i = ipivot(k)
         ipivot(k) = ipivot(jbar)
         ipivot(jbar) = i
         sum(jbar) = sum(k)
         sum(k) = sigma
c
         do 40 i=1,m
              sigma = qr(i,k)
              qr(i,k) = qr(i,jbar)
              qr(i,jbar) = sigma
   40    continue
c        ..........end of column interchange..........
   50    continue
c        ..........  second inner product  ..........
         qrkmax = zero
c
         do 60 i=k,m
              if (abs( qr(i,k) ) .gt. qrkmax)  qrkmax = abs( qr(i,k) )
   60    continue
c
         if (qrkmax .lt. ufeta)  go to 210
         alphak = v2norm(mk1, qr(k,k)) / qrkmax
         sigma = alphak**2
c
c        ..........  end second inner product  ..........
         qrkk = qr(k,k)
         if (qrkk .ge. zero)  alphak = -alphak
         alpha(k) = alphak * qrkmax
         beta = qrkmax * sqrt(sigma - (qrkk*alphak/qrkmax) )
         qr(k,k) = qrkk - alpha(k)
         do 65 i=k,m
   65         qr(i,k) =  qr(i,k) / beta
         k1 = k + 1
         if (k1 .gt. n) go to 120
c
         do 110 j = k1, n
              temp = -dotprd(mk1, qr(k,k), qr(k,j))
c
c             ***  set qr(i,j) = qr(i,j) + temp*qr(i,k), i = k,...,m.
c
              call vaxpy(mk1, qr(k,j), temp, qr(k,k), qr(k,j))
c
              if (k1 .gt. m) go to 110
              sumj = sum(j)
              if (sumj .lt. ufeta) go to 110
              temp = abs(qr(k,j)/sumj)
              if (temp .lt. rktol1) go to 110
              if (temp .ge. p99) go to 90
                   sum(j) = sumj * sqrt(one - temp**2)
                   go to 110
   90         sum(j) = v2norm(m-k, qr(k1,j))
  110    continue
c        ..........end of k-th householder transformation..........
  120 continue
c
      go to 999
c     ..........error exit on k-th transformation..........
  210 ierr = -k
      go to 230
c     ..........no non-zero acceptable pivot found..........
  220 ierr = k
  230 do 240 i = k, n
         alpha(i) = zero
         if (i .gt. k) call vscopy(i-k, qr(k,i), zero)
 240     continue
c     ..........return to caller..........
  999 return
c     ..........last card of qrfact..........
      end
      real function reldst(p, d, x, x0)
c
c  ***  compute and return relative difference between x and x0  ***
c  ***  nl2sol version 2.2  ***
c
      integer p
      real d(p), x(p), x0(p)
c/+
      real abs
c/
      integer i
      real emax, t, xmax, zero
c/6
      data zero/0.e+0/
c/7
c     parameter (zero=0.d+0)
c/
c
      emax = zero
      xmax = zero
      do 10 i = 1, p
         t = abs(d(i) * (x(i) - x0(i)))
         if (emax .lt. t) emax = t
         t = d(i) * (abs(x(i)) + abs(x0(i)))
         if (xmax .lt. t) xmax = t
 10      continue
      reldst = zero
      if (xmax .gt. zero) reldst = emax / xmax
 999  return
c  ***  last card of reldst follows  ***
      end
      subroutine rptmul(func, ipivot, j, nn, p, rd, x, y, z)
c
c  ***  func = 1... set  y = rmat * (perm**t) * x.
c  ***  func = 2... set  y = perm * (rmat**t) * rmat * (perm**t) * x.
c  ***  func = 3... set  y = perm * (rmat**t) x.
c
c
c  ***  perm = matrix whose i-th col. is the ipivot(i)-th unit vector.
c  ***  rmat is the upper triangular matrix whose strict upper triangle
c  ***       is stored in  j  and whose diagonal is stored in rd.
c  ***  z is a scratch vector.
c  ***  x and y may share storage.
c
      integer func, nn, p
      integer ipivot(p)
      real j(nn,p), rd(p), x(p), y(p), z(p)
c
c  ***  local variables  ***
c
      integer i, im1, k, km1
      real zk
c
c  ***  external function  ***
c
      external dotprd
      real dotprd
c
c-----------------------------------------------------------------------
c
      if (func .gt. 2) go to 50
c
c  ***  first set  z = (perm**t) * x  ***
c
      do 10 i = 1, p
         k = ipivot(i)
         z(i) = x(k)
 10      continue
c
c  ***  now set  y = rmat * z  ***
c
      y(1) = z(1) * rd(1)
      if (p .le. 1) go to 40
      do 30 k = 2, p
         km1 = k - 1
         zk = z(k)
         do 20 i = 1, km1
 20           y(i) = y(i) + j(i,k)*zk
         y(k) = zk*rd(k)
 30      continue
c
 40   if (func .le. 1) go to 999
      go to 70
c
 50   do 60 i = 1, p
 60      y(i) = x(i)
c
c  ***  set  z = (rmat**t) * y  ***
c
 70   z(1) = y(1) * rd(1)
      if (p .eq. 1) go to 90
      do 80 i = 2, p
         im1 = i - 1
         z(i) = y(i)*rd(i) + dotprd(im1, j(1,i), y)
 80      continue
c
c  ***  now set  y = perm * z  ***
c
 90   do 100 i = 1, p
         k = ipivot(i)
         y(k) = z(i)
 100     continue
c
 999  return
c  ***  last card of rptmul follows  ***
      end
      subroutine slupdt(a, cosmin, p, size, step, u, w, wchmtd, wscale,
     1                  y)
c
c  ***  update symmetric  a  so that  a * step = y  ***
c  ***  (lower triangle of  a  stored rowwise       ***
c
c  ***  parameter declarations  ***
c
      integer p
      real a(1), cosmin, size, step(p), u(p), w(p),
     1                 wchmtd(p), wscale, y(p)
c     dimension a(p*(p+1)/2)
c
c  ***  local variables  ***
c
      integer i, j, k
      real denmin, sdotwm, t, ui, wi
c
c     ***  constants  ***
      real half, one, zero
c
c  ***  intrinsic functions  ***
c/+
      real abs, amin1
c/
c  ***  external functions and subroutines  ***
c
      external dotprd, slvmul, v2norm
      real dotprd, v2norm
c
c/6
      data half/0.5e+0/, one/1.e+0/, zero/0.e+0/
c/7
c     parameter (half=0.5d+0, one=1.d+0, zero=0.d+0)
c/
c
c-----------------------------------------------------------------------
c
      sdotwm = dotprd(p, step, wchmtd)
      denmin = cosmin * v2norm(p,step) * v2norm(p,wchmtd)
      wscale = one
      if (denmin .ne. zero) wscale = amin1(one, abs(sdotwm/denmin))
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
      subroutine slvmul(p, y, s, x)
c
c  ***  set  y = s * x,  s = p x p symmetric matrix.  ***
c  ***  lower triangle of  s  stored rowwise.         ***
c
c  ***  parameter declarations  ***
c
      integer p
      real s(1), x(p), y(p)
c     dimension s(p*(p+1)/2)
c
c  ***  local variables  ***
c
      integer i, im1, j, k
      real xi
c
c  ***  no intrinsic functions  ***
c
c  ***  external function  ***
c
      external dotprd
      real dotprd
c
c-----------------------------------------------------------------------
c
      j = 1
      do 10 i = 1, p
         y(i) = dotprd(i, s(j), x)
         j = j + i
 10      continue
c
      if (p .le. 1) go to 999
      j = 1
      do 40 i = 2, p
         xi = x(i)
         im1 = i - 1
         j = j + 1
         do 30 k = 1, im1
              y(k) = y(k) + s(j)*xi
              j = j + 1
 30           continue
 40      continue
c
 999  return
c  ***  last card of slvmul follows  ***
      end
      logical function stopx(idummy)
c     *****parameters...
      integer idummy
c
c     ..................................................................
c
c     *****purpose...
c     this function may serve as the stopx (asynchronous interruption)
c     function for the nl2sol (nonlinear least-squares) package at
c     those installations which do not wish to implement a
c     dynamic stopx.
c
c     *****algorithm notes...
c     at installations where the nl2sol system is used
c     interactively, this dummy stopx should be replaced by a
c     function that returns .true. if and only if the interrupt
c     (break) key has been pressed since the last call on stopx.
c
c     ..................................................................
c
      stopx = .false.
      return
      end
      subroutine vaxpy(p, w, a, x, y)
c
c  ***  set w = a*x + y  --  w, x, y = p-vectors, a = scalar  ***
c
      integer p
      real a, w(p), x(p), y(p)
c
      integer i
c
      do 10 i = 1, p
 10      w(i) = a*x(i) + y(i)
      return
      end
      subroutine vcopy(p, y, x)
c
c  ***  set y = x, where x and y are p-vectors  ***
c
      integer p
      real x(p), y(p)
c
      integer i
c
      do 10 i = 1, p
 10      y(i) = x(i)
      return
      end
      subroutine vscopy(p, y, s)
c
c  ***  set p-vector y to scalar s  ***
c
      integer p
      real s, y(p)
c
      integer i
c
      do 10 i = 1, p
 10      y(i) = s
      return
      end
      real function v2norm(p, x)
c
c  ***  return the 2-norm of the p-vector x, taking  ***
c  ***  care to avoid the most likely underflows.    ***
c
      integer p
      real x(p)
c
      integer i, j
      real one, r, scale, sqteta, t, xi, zero
c/+
      real abs, sqrt
c/
      external rmdcon
      real rmdcon
c
c/6
      data one/1.e+0/, zero/0.e+0/
c/7
c     parameter (one=1.d+0, zero=0.d+0)
c     save sqteta
c/
      data sqteta/0.e+0/
c
      if (p .gt. 0) go to 10
         v2norm = zero
         go to 999
 10   do 20 i = 1, p
         if (x(i) .ne. zero) go to 30
 20      continue
      v2norm = zero
      go to 999
c
 30   scale = abs(x(i))
      if (i .lt. p) go to 40
         v2norm = scale
         go to 999
 40   t = one
      if (sqteta .eq. zero) sqteta = rmdcon(2)
c
c     ***  sqteta is (slightly larger than) the square root of the
c     ***  smallest positive floating point number on the machine.
c     ***  the tests involving sqteta are done to prevent underflows.
c
      j = i + 1
      do 60 i = j, p
         xi = abs(x(i))
         if (xi .gt. scale) go to 50
              r = xi / scale
              if (r .gt. sqteta) t = t + r*r
              go to 60
 50           r = scale / xi
              if (r .le. sqteta) r = zero
              t = one  +  t * r*r
         scale = xi
 60      continue
c
      v2norm = scale * sqrt(t)
 999  return
c  ***  last card of v2norm follows  ***
      end
c///////////////////////////////////////////////////////////////////////
c  ***  run nl2sol on various test problems, print summary statistics.
c
c     *****common storage with nltest.
c
      common /testcm/ v, rs, jac, nout, nprob, xscal1, xscal2, is, iv
      common /testch/ name, irc
      integer is(6,50), iv(80), jac, nout, nprob, xscal1, xscal2
      real rs(5,50)
c/6
      real name(2,50)
      integer irc(50)
c/7
c     character name(2,50)*4, irc(50)*1
c/
      real v(1736)
c
c
c     ..................................................................
c
c     *****purpose.
c        this main program calls nltest to run nl2sol, the nonlinear
c     least-squares solver of ref. 1, on various test problems.
c
c
c     *****application and usage restrictions.
c     this main driver is intended to check whether the nl2sol
c     (nonlinear least-squares) package was successfully
c     transported to a new machine.
c
c     *****algorithm notes.
c     the test problems used are from references (2), (3), and (4).
c     some additional test problems were suggested by jorge more (pri-
c     vate communication).  calls passing these problems to nltest have
c     been commented out (since there are enough other problems), but
c     not removed, since they may be of interest to other researchers.
c
c     *****functions and subroutines called.
c
c        dfault - establishes the default parameter settings for
c                 iv and v.
c
c        imdcon - imdcon(2) returns i/o unit number on which nltest
c                  writes a summary of each test run.
c
c        ivvset - supplies nondefault values for iv and v.
c
c        nltest - calls nl2sol, the nonlinear least-squares
c                  problem solver.
c
c        today  - supplies date and time (or current version of nl2sol).
c
c     *****references.
c
c     (1). dennis, j.e.. gay, d.m.. and welsch, r.e. (1980),
c          an adaptive nonlinear least-squares algorithm,
c          submitted to acm trans. math. software.
c          under revision.
c
c     (2). gill, p.e.. and murray, w. (1976),algorithms for the
c          solution of the non-linear least-squares problem,
c          npl report nac71,(national physical laboratory,
c          division of numerical analysis and computing,
c          teddington,middlesex,england).
c
c     (3) meyer, r.r. (1970), theoretical and computational aspects
c        of nonlinear regression, pp. 465-486 of nonlinear programming,
c        edited by j.b. rosen, o.l.mangasarian, and k. ritter,
c        academic press, new york.
c
c     (4) brown, k.m. (1969), a quadratically convergent newton-
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
c     *****intrinsic functions.
c/+
      integer mod
      real amax1
c/
c     *****external functions and subroutines.
      external dfault, imdcon, ivvset, nltest, today
      integer imdcon
c
c     *****local variables.
      logical rstart
      integer i, j, k, mxfcsv, mxitsv, pu
c/6
      integer jtyp(2)
      real datime(4)
c/7
c     character datime(4)*4, jtyp(2)*1
c/
c
c/6
      data rstart/.false./, jtyp(1),jtyp(2)/1h ,1h*/
c/7
c     data rstart/.false./, jtyp(1),jtyp(2)/' ','*'/
c/
c
c-----------------------------------------------------------------------
c
c  ***  establish default parameter settings  ***
      call dfault (iv, v)
      nout = imdcon(2)
c
c  ***  non-default parameter settings  ***
c
      call ivvset(iv, v)
      pu = iv(21)
c
      jac = 1
      nprob = 0
      xscal1 = 1
      xscal2 = 3
c
c/6
      call nltest(2,2,1,4hrosn,4hbrok,rstart)
      call nltest(3,3,2,4hheli,4hx   ,rstart)
      call nltest(4,4,3,4hsing,4hular,rstart)
      call nltest(7,4,4,4hwood,4hs   ,rstart)
      xscal2 = 1
      call nltest(3,3,5,4hzang,4hwill,rstart)
      xscal2 = 3
      call nltest(5,3,6,4hengv,4hall ,rstart)
      call nltest(2,2,7,4hbran,4hin  ,rstart)
      xscal2 = 2
      call nltest(3,2,8,4hbeal,4he   ,rstart)
      call nltest(5,4,9,4hcrag,4hg   ,rstart)
      xscal2 = 2
      call nltest(10,3,10,4hbox ,4h    ,rstart)
      mxfcsv = iv(17)
      mxitsv = iv(18)
      iv(17) = 20
      iv(18) = 15
      xscal2 = 1
      call nltest(15,15,11,4hdavi,4hdon1,rstart)
      iv(17) = mxfcsv
      iv(18) = mxitsv
      xscal2 = 3
      call nltest(2,2,12,4hfrds,4htein,rstart)
      xscal2 = 1
      call nltest(31,6,13,4hwats,4hon6 ,rstart)
      call nltest(31,9,14,4hwats,4hon9 ,rstart)
      call nltest(31,12,15,4hwats,4hon12,rstart)
      mxfcsv = iv(17)
      iv(17) = 20
      mxitsv = iv(18)
      iv(18) = 15
      call nltest(31,20,16,4hwats,4hon20,rstart)
      iv(17) = mxfcsv
      iv(18) = mxitsv
      xscal2 = 2
      call nltest(8,8,17,4hcheb,4hqd8 ,rstart)
      xscal2 = 3
      call nltest(20,4,18,4hbrow,4hn   ,rstart)
      call nltest(15,3,19,4hbard,4h    ,rstart)
      xscal2 = 1
      call nltest(10,2,20,4hjenn,4hrich,rstart)
      xscal2 = 3
      call nltest(11,4,21,4hkowa,4hlik ,rstart)
      xscal2 = 1
      call nltest(33,5,22,4hosbo,4hrne1,rstart)
      xscal2 = 2
      call nltest(65,11,23,4hosbo,4hrne2,rstart)
      xscal2 = 3
      call nltest(3,2,24,4hmads,4hen  ,rstart)
      xscal2 = 1
      iv(17) = 400
      iv(18) = 300
      call nltest(16,3,25,4hmeye,4hr   ,rstart)
c/7
c     call nltest(2,2,1,'rosn','brok',rstart)
c     call nltest(3,3,2,'heli','x   ',rstart)
c     call nltest(4,4,3,'sing','ular',rstart)
c     call nltest(7,4,4,'wood','s   ',rstart)
c     xscal2 = 1
c     call nltest(3,3,5,'zang','will',rstart)
c     xscal2 = 3
c     call nltest(5,3,6,'engv','all ',rstart)
c     call nltest(2,2,7,'bran','in  ',rstart)
c     xscal2 = 2
c     call nltest(3,2,8,'beal','e   ',rstart)
c     call nltest(5,4,9,'crag','g   ',rstart)
c     xscal2 = 2
c     call nltest(10,3,10,'box ','    ',rstart)
c     mxfcsv = iv(17)
c     mxitsv = iv(18)
c     iv(17) = 20
c     iv(18) = 15
c     xscal2 = 1
c     call nltest(15,15,11,'davi','don1',rstart)
c     iv(17) = mxfcsv
c     iv(18) = mxitsv
c     xscal2 = 3
c     call nltest(2,2,12,'frds','tein',rstart)
c     xscal2 = 1
c     call nltest(31,6,13,'wats','on6 ',rstart)
c     call nltest(31,9,14,'wats','on9 ',rstart)
c     call nltest(31,12,15,'wats','on12',rstart)
c     mxfcsv = iv(17)
c     iv(17) = 20
c     mxitsv = iv(18)
c     iv(18) = 15
c     call nltest(31,20,16,'wats','on20',rstart)
c     iv(17) = mxfcsv
c     iv(18) = mxitsv
c     xscal2 = 2
c     call nltest(8,8,17,'cheb','qd8 ',rstart)
c     xscal2 = 3
c     call nltest(20,4,18,'brow','n   ',rstart)
c     call nltest(15,3,19,'bard','    ',rstart)
c     xscal2 = 1
c     call nltest(10,2,20,'jenn','rich',rstart)
c     xscal2 = 3
c     call nltest(11,4,21,'kowa','lik ',rstart)
c     xscal2 = 1
c     call nltest(33,5,22,'osbo','rne1',rstart)
c     xscal2 = 2
c     call nltest(65,11,23,'osbo','rne2',rstart)
c     xscal2 = 3
c     call nltest(3,2,24,'mads','en  ',rstart)
c     xscal2 = 1
c     iv(17) = 400
c     iv(18) = 300
c     call nltest(16,3,25,'meye','r   ',rstart)
c/
c  ***  brown5  ***
c     call nltest(5,5,26,4hbrow,4hn5  ,rstart)
c  ***  brown10  ***
c     call nltest(10,10,27,4hbrow,4hn10 ,rstart)
c  ***  brown30  ***
c     call nltest(30,30,28,4hbrow,4hn30 ,rstart)
c  ***  brown40  ***
c     call nltest(40,40,29,4hbrow,4hn40 ,rstart)
c  ***  bard+10 ***
c     call nltest(15,3,30,4hbard,4h+10 ,rstart)
c  ***  kowalik and osborne + 10  ***
c     call nltest(11,4,31,4hkowa,4hl+10,rstart)
c  ***  meyer + 10  ***
c     call nltest(16,3,32,4hmeye,4hr+10,rstart)
c  ***  watson6 + 10  ***
c     call nltest(31,6,33,4hwat6,4h+10 ,rstart)
c  ***  watson9 + 10  ***
c     call nltest(31,9,34,4hwat9,4h+10 ,rstart)
c  ***  watson12 + 10  ***
c     call nltest(31,12,35,4hwat1,4h2+10,rstart)
c  ***  watson20 + 10  ***
c     call nltest(31,20,36,4hwat2,4h0+10,rstart)
c
c  ***  repeat two tests using finite-difference jacobian  ***
c
      jac = 2
      xscal2 = 1
c
      iv(17) = 50
      iv(18) = 40
c/6
      call nltest(2,2,1,4hrosn,4hbrok,rstart)
c/7
c     call nltest(2,2,1,'rosn','brok',rstart)
c/
      v(29) = amax1(1.0e-7, v(29))
      iv(17) = 30
      iv(18) = 20
c  ***  brown  ***
c/6
      call nltest(20,4,18,4hbrow,4hn   ,rstart)
c/7
c     call nltest(20,4,18,'brow','n   ',rstart)
c/
c
      if (nprob .eq. 0 .or. pu .eq. 0) stop
      call today(datime)
      do 130 k = 1, nprob
         if (mod(k,56) .eq. 1) write(pu, 110) datime, nprob
 110     format(1h1,11x,2a4,2x,2a4,10x,10hsummary of,i4,
     1          22h nl2sol test runs.....,10x,
     2          32h(* = finite-difference jacobian)/
     3          48h0 problem    n   p  niter   nf   ng  iv1  x0scal,5x,
     4          39hfinal f     preldf     nreldf     reldx/)
         j = is(6,k)
         write(pu,120) jtyp(j), name(1,k), name(2,k),
     1                 (is(i,k), i=1,5), irc(k), (rs(i,k), i=1,5)
 120     format(1x,a1,2a4,2i4,i7,2i5,3x,a1,f9.1,e13.3,3e11.3)
 130     continue
c
      stop
c...... last card of nlmain ............................................
      end
      subroutine ivvset(iv, v)
c
c  ***  supply nondefault iv and v values for nlmain  (nl2sol ver. 2.2).
c
      integer iv(24)
      real v(100)
c
c     activate the next line to turn off detailed summary printing
c     iv(21) = 0
      return
      end
      subroutine nltest (n, p, nex, title1, title2, rstart)
c
c  ***  call nl2sol, save and print statistics  ***
c
c
      integer n, p, nex
      logical rstart
c/6
      real title1, title2
c/7
c     character*4 title1, title2
c/
c
      common /testcm/ v, rs, jac, nout, nprob, xscal1, xscal2, is, iv
      common /testch/ name, irc
      integer is(6,50), iv(80), jac, nout, nprob, xscal1, xscal2
      real rs(5,50)
c/6
      integer irc(50)
      real name(2,50)
c/7
c     character name(2,50)*4, irc(50)*1
c/
      real v(1736)
c
      logical rstrt
      integer i, irun, pu, uip(1)
c/6
      integer alg(2), jtyp(2), rc(10)
      real datime(4)
c/7
c     character*4 datime(4)
c     character*2 alg(2)
c     character*1 jtyp(2), rc(10)
c/
      real one, t, urparm(1), x(20), x0scal, zero
c
c     ***  external functions and subroutines  ***
c
      external nl2sno, nl2sol, testr, testj, today, xinit
c
c  ***  iv and v subscripts  ***
c
      integer f, f0, nfcall, nfcov, ngcall, niter, nreduc, preduc,
     1        prunit, reldx
c
c/6
      data f/10/, f0/13/, nfcall/6/, nfcov/40/, ngcall/30/,
     1     ngcov/41/, niter/31/, nreduc/6/, preduc/7/,
     2     prunit/21/, reldx/17/
c/7
c     parameter (f=10, f0=13, nfcall=6, nfcov=40, ngcall=30,
c    1     ngcov=41, niter=31, nreduc=6, preduc=7,
c    2     prunit=21, reldx=17)
c/
c/6
      data one/1.e+0/, zero/0.e+0/
c/7
c     parameter (one=1.d+0, zero=0.d+0)
c/
c/6
      data alg(1),alg(2)/2hol,2hno/, jtyp(1),jtyp(2)/1h ,1h*/
      data rc(1)/1h./, rc(2)/1h+/, rc(3)/1hx/, rc(4)/1hr/, rc(5)/1hb/,
     1     rc(6)/1ha/, rc(7)/1hs/, rc(8)/1hf/, rc(9)/1he/, rc(10)/1hi/
c/7
c     data alg(1),alg(2)/'ol','no'/, jtyp(1),jtyp(2)/' ','*'/
c     data rc(1)/'.'/, rc(2)/'+'/, rc(3)/'x'/, rc(4)/'r'/, rc(5)/'b'/,
c    1     rc(6)/'a'/, rc(7)/'s'/, rc(8)/'f'/, rc(9)/'e'/, rc(10)/'i'/
c/
c
c-----------------------------------------------------------------------
c
      uip(1) = nex
      rstrt = rstart
      if (rstrt) go to 20
         pu = iv(prunit)
         call today(datime)
         if (pu .ne. 0) write(pu,10) alg(jac), title1, title2, datime
 10      format (1h1//11h ***** nl2s,a2,12h on problem ,2a4,6h *****,6x,
     1           2a4,2x,2a4)
c
 20   do 100 irun = xscal1, xscal2
         if (rstrt) go to 40
         iv(1) = 12
         x0scal = 1.0e1 ** (irun-1)
c
c        ***  initialize the solution vector x  ***
         call xinit(p, x, nex)
         do 30 i = 1, p
 30           x(i) = x0scal * x(i)
c
 40      if (jac .eq. 1)
     1             call nl2sol(n,p,x,testr,testj,iv,v,uip,urparm,testr)
         if (jac .eq. 2)
     1             call nl2sno(n,p,x,testr,iv,v,uip,urparm,testr)
         if (.not. rstrt .and. nprob .lt. 50) nprob = nprob + 1
         name(1,nprob) = title1
         name(2,nprob) = title2
         is(1,nprob) = n
         is(2,nprob) = p
         is(3,nprob) = iv(niter)
         is(4,nprob) = iv(nfcall) - iv(nfcov)
         is(5,nprob) = iv(ngcall) - iv(ngcov)
         i = iv(1)
         irc(nprob) = rc(i)
         is(6,nprob) = jac
         rs(1,nprob) = x0scal
         rs(2,nprob) = v(f)
         t = one
         if (v(f0) .gt. zero) t = v(preduc) / v(f0)
         rs(3,nprob) = t
         t = one
         if (v(f0) .gt. zero) t = v(nreduc) / v(f0)
         rs(4,nprob) = t
         rs(5,nprob) = v(reldx)
         rstrt = .false.
         if (nout .eq. 0) go to 100
         if (nprob .eq. 1) write(nout,50) datime
 50      format(1h1,11x,2a4,2x,2a4,10x,24hnl2sol test summary.....,10x,
     1          32h(* = finite-difference jacobian)/
     2          48h0 problem    n   p  niter   nf   ng  iv1  x0scal,5x,
     3          39hfinal f     preldf     nreldf     reldx/)
         write(nout,60) jtyp(jac), title1, title2,
     1                (is(i,nprob),i=1,5),irc(nprob),(rs(i,nprob),i=1,5)
 60      format(1x,a1,2a4,2i4,i7,2i5,3x,a1,f9.1,e13.3,3e11.3)
 100     continue
c
 999  return
c  ***  last card of nltest follows  ***
      end
      subroutine testj(n, p, x, nfcall, j, uiparm, urparm, ufparm)
c
c  ***  parameters  ***
c
      integer n, p, nfcall, uiparm(1)
      real x(p), j(n,p), urparm(1)
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
      real e, expmin, r2, t, theta, ti, tim1, tip1, tpi,
     1   tpim1, tpip1, twopi, u, uftolg, ukow(11), v, w, z, zero
      integer i, k, nex, nm1
c  ***  intrinsic functions  ***
c/+
      real alog, amin1, cos, exp, float, sin, sqrt
c/
      external rmdcon
      real rmdcon
c
c/6
c /6
      data twopi/6.283185e+0/, zero/0.e+0/
c /7
c     parameter (twopi=6.283185e+0, zero=0.e+0)
c /
c/6
c/7
c     save expmin, uftolg
c/
      data ukow(1)/4.0/, ukow(2)/2.0/, ukow(3)/1.0/,
     1   ukow(4)/5.0e-1/, ukow(5)/2.5e-1/, ukow(6)/1.67e-1/,
     2   ukow(7)/1.25e-1/, ukow(8)/1.0e-1/, ukow(9)/8.33e-2/,
     3   ukow(10)/7.14e-2/, ukow(11)/6.25e-2/
c  ***  machine dependent constant  ***
      data expmin/0.0/, uftolg/0./
c
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
 100  j(1,1) = -2.0e1*x(1)
      j(1,2) = 1.0e1
      j(2,1) = -1.0
      j(2,2) = 0.0
      go to 9999
c  ***  helix  ***
 200  t = x(1)**2 + x(2)**2
      ti = 1.e2/(twopi*t)
      j(1,1) = ti*x(2)
      t = 1.e1/sqrt(t)
      j(2,1) = x(1)*t
      j(3,1) = 0.
      j(1,2) = -ti*x(1)
      j(2,2) = x(2)*t
      j(3,2) = 0.
      j(1,3) = 1.e1
      j(2,3) = 0.
      j(3,3) = 1.
      go to 9999
c  ***  singular  ***
 300  do 301 k = 1,4
         do 301 i = 1,4
 301          j(i,k) = 0.
      j(1,1) = 1.
      j(1,2) = 1.e1
      j(2,3) = sqrt(5.)
      j(2,4) = -j(2,3)
      j(3,2) = 2.*(x(2) - 2.*x(3))
      j(3,3) = -2.*j(3,2)
      j(4,1) = sqrt(4.e1)*(x(1) - x(4))
      j(4,4) = -j(4,1)
      go to 9999
c  ***  woods  ***
 400  do 401 k = 1,4
         do 401 i = 1,7
 401            j(i,k) = 0.
      j(1,1) = -2.e1*x(1)
      j(1,2) = 1.e1
      j(2,1) = -1.
      j(3,4) = sqrt(9.e1)
      j(3,3) = -2.*x(3)*j(3,4)
      j(4,3) = -1.
      j(5,2) = sqrt(9.9)
      j(5,4) = j(5,2)
      j(6,2) = sqrt(0.2)
      j(7,4) = j(6,2)
      go to 9999
c  ***  zangwill  ***
 500  do 501 k = 1,3
         do 501 i = 1,3
 501            j(i,k) = 1.
      j(1,2) = -1.
      j(2,1) = -1.
      j(3,3) = -1.
      go to 9999
c  ***  engvall  ***
 600  j(1,1) = 2.*x(1)
      j(1,2) = 2.*x(2)
      j(1,3) = 2.*x(3)
      j(2,1) = j(1,1)
      j(2,2) = j(1,2)
      j(2,3) = 2.*(x(3) - 2.)
      j(3,1) = 1.
      j(3,2) = 1.
      j(3,3) = 1.
      j(4,1) = 1.
      j(4,2) = 1.
      j(4,3) = -1.
      t = 2.*(5.*x(3) - x(1) + 1.)
      j(5,1) = 3.*x(1)**2 - t
      j(5,2) = 6.*x(2)
      j(5,3) = 5.*t
      go to 9999
c  ***  branin  ***
 700  j(1,1) = 4.
      j(1,2) = 4.
      j(2,1) = 3. + (x(1) - 2.)*(3.*x(1) - 2.*x(2) - 2.) +
     1   x(2)*x(2)
      j(2,2) = 1. + 2.*(2.*x(1) - x(2)*x(2)) - (x(1) - x(2))**2
      go to 9999
c  ***  beale  ***
 800  j(1,1) = x(2) - 1.
      j(1,2) = x(1)
      j(2,1) = x(2)**2 - 1.
      j(2,2) = 2.*x(1)*x(2)
      j(3,1) = x(2)**3 - 1.
      j(3,2) = 3.*x(1)*(x(2)**2)
      go to 9999
c  ***  cragg & levy  ***
 900  do 901 i = 1,5
         do 901 k = 1,4
 901          j(i,k) = 0.
      t = exp(x(1))
      j(1,2) = -2.*(t - x(2))
      j(1,1) = -t * j(1,2)
      j(2,2) = 3.0e1*(x(2) - x(3))**2
      j(2,3) = -j(2,2)
      j(3,3) = 2.*sin(x(3) - x(4))/(cos(x(3) - x(4)))**3
      j(3,4) = -j(3,3)
      j(4,1) = 4.*x(1)**3
      j(5,4) = 1.
      go to 9999
c  ***  box  ***
 1000 if (expmin .eq. zero) expmin = 1.999*alog(rmdcon(2))
      do 1001 i = 1,10
         ti = -0.1*float(i)
         e = zero
         t = x(1)*ti
         if (t .ge. expmin) e = exp(t)
         j(i,1) = ti*e
         e = zero
         t = x(2)*ti
         if (t .ge. expmin) e = exp(t)
         j(i,2) = -ti*e
         j(i,3) = exp(1.e1*ti) - exp(ti)
 1001    continue
      go to 9999
c  ***  davidon 1  ***
 1100 nm1 = n-1
      do 1101 i = 1,nm1
         ti = float(i)
         t = 1.
         do 1101 k = 1,p
              j(i,k) = t
              t = t*ti
 1101         continue
      j(n,1) = 1.
      do 1102 k = 2,p
 1102    j(n,k) = 0.
      go to 9999
c  ***  freudenstein & roth  ***
 1200 j(1,1) = 1.
      j(1,2) = -2. + x(2)*(1.e1 - 3.*x(2))
      j(2,1) = 1.
      j(2,2) = -1.4e1 + x(2)*(2. + 3.*x(2))
      go to 9999
c  ***  watson  ***
 1300 continue
 1400 continue
 1500 continue
 1600 do 1603 i = 1,29
         ti = float(i)/2.9e1
         r2 = x(1)
         t= 1.
         do 1601 k = 2,p
              t = t*ti
              r2 = r2 + t*x(k)
 1601    continue
         r2 = -2.*r2
         j(i,1) = r2
         t = 1.
         r2 = ti*r2
         do 1602 k = 2,p
              j(i,k) = t*(float(k-1) + r2)
              t = t*ti
 1602    continue
 1603 continue
      do 1604 i = 30,31
         do 1604 k = 2,p
 1604         j(i,k) = 0.
      j(30,1) = 1.
      j(31,1) = -2.*x(1)
      j(31,2) = 1.
      go to 9999
c  ***  chebyquad  ***
 1700 do 1701 k = 1,n
         tim1 = -1./float(n)
         z = 2.*x(k) - 1.
         ti = z*tim1
         tpim1 = 0.
         tpi = 2.*tim1
         z = z + z
         do 1701 i = 1,n
              j(i,k) = tpi
              tpip1 = 4.*ti + z*tpi - tpim1
              tpim1 = tpi
              tpi = tpip1
              tip1 = z*ti - tim1
              tim1 = ti
              ti = tip1
 1701         continue
      go to 9999
c  ***  brown and dennis  ***
 1800 do 1801 i = 1, n
         ti = 0.2*float(i)
         j(i,1) = 2.0*(x(1) + x(2)*ti - exp(ti))
         j(i,2) = ti*j(i,1)
         t = sin(ti)
         j(i,3) = 2.0*(x(3) + x(4)*t - cos(ti))
         j(i,4) = t*j(i,3)
 1801    continue
      go to 9999
c  ***  bard  ***
 1900 do 1901 i = 1,15
         j(i,1) = -1.
         u = float(i)
         v = 1.6e1 - u
         w = amin1 (u,v)
         t = u/(x(2)*v + x(3)*w)**2
         j(i,2) = v*t
         j(i,3) = w*t
 1901 continue
      go to 9999
c  *** jennrich & sampson  ***
 2000 do 2001 i = 1,10
         ti = float(i)
         j(i,1) = -ti*exp(ti*x(1))
         j(i,2) = -ti*exp(ti*x(2))
 2001    continue
      go to 9999
c  ***  kowalik & osborne  ***
 2100 do 2101 i = 1,11
         t = -1./(ukow(i)**2 + x(3)*ukow(i) + x(4))
         j(i,1) = t*(ukow(i)**2 + x(2)*ukow(i))
         j(i,2) = x(1)*ukow(i)*t
         t = t*j(i,1)*x(1)
         j(i,3) = ukow(i)*t
         j(i,4) = t
 2101 continue
      go to 9999
c  ***  osborne 1  ***
 2200 do 2201 i = 1,33
         ti = 1.0e1*float(1-i)
         j(i,1) = -1.
         j(i,2) = -exp(x(4)*ti)
         j(i,3) = -exp(x(5)*ti)
         j(i,4) = ti*x(2)*j(i,2)
         j(i,5) = ti*x(3)*j(i,3)
 2201    continue
      go to 9999
c  ***  osborne 2  ***
c     ***  uftolg is a machine-dependent constant.  it is just slightly
c     ***  larger than the log of the smallest positive machine number.
 2300 if (uftolg .eq. 0.) uftolg = 1.999 * alog(rmdcon(2))
      do 2302 i = 1,65
         ti = float(1 - i)*1.e-1
         j(i,1) = -exp(x(5)*ti)
         j(i,5) = x(1)*ti*j(i,1)
         do 2301 k = 2,4
              t = x(k + 7) + ti
              r2 = 0.
              theta = -x(k+4)*t*t
              if (theta .gt. uftolg) r2 = -exp(theta)
              j(i,k) = r2
              r2 = -t*r2*x(k)
              j(i,k+4) = r2*t
              j(i,k+7) = 2.*x(k+4)*r2
 2301         continue
 2302    continue
      go to 9999
c  ***  madsen  ***
 2400 j(1,1) = 2.*x(1) + x(2)
      j(1,2) = 2.*x(2) + x(1)
      j(2,1) = cos(x(1))
      j(2,2) = 0.
      j(3,1) = 0.
      j(3,2) = -sin(x(2))
      go to 9999
c  ***  meyer  ***
 2500 do 2501 i = 1, 16
         ti = float(5*i + 45)
         u = ti + x(3)
         t = exp(x(2)/u)
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
              j(i,k) = 1.0
              if (i .eq. k) j(i,k) = 2.0
 2901         continue
      do 2903 k = 1, n
         t = 1.0
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
      subroutine testr(n, p, x, nfcall, r, uiparm, urparm, ufparm)
c
c     *****parameters.
c
      integer n, p, nfcall, uiparm(1)
      real x(p), r(n), urparm(1)
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
      real e1, e2, floatn, ri, r1, r2, t, theta, ti, tim1,
     1             tip1, twopi, t1, t2, u, v, w, z
      real ybard(15), ykow(11), ukow(11), yosb1(33),
     1             yosb2(65), ymeyer(16)
      integer i, j, nex, nm1
      real expmax, expmin, uftolg
c  ***  intrinsic functions  ***
c/+
      integer mod
      real alog, amin1, atan2, cos, exp, float, sin, sqrt
c/
      external rmdcon
      real rmdcon
c /6
      data twopi/6.283185e+0/
c /7
c     parameter (twopi=6.283185e+0)
c /
c/6
c/7
c     save expmax, expmin, uftolg
c/
      data ybard(1)/1.4e-1/, ybard(2)/1.8e-1/, ybard(3)/2.2e-1/,
     1   ybard(4)/2.5e-1/, ybard(5)/2.9e-1/, ybard(6)/3.2e-1/,
     2   ybard(7)/3.5e-1/, ybard(8)/3.9e-1/, ybard(9)/3.7e-1/,
     3   ybard(10)/5.8e-1/, ybard(11)/7.3e-1/, ybard(12)/9.6e-1/,
     4   ybard(13)/1.34/, ybard(14)/2.10/, ybard(15)/4.39/
      data ykow(1)/1.957e-1/, ykow(2)/1.947e-1/, ykow(3)/1.735e-1/,
     1   ykow(4)/1.600e-1/, ykow(5)/8.44e-2/, ykow(6)/6.27e-2/,
     2   ykow(7)/4.56e-2/, ykow(8)/3.42e-2/, ykow(9)/3.23e-2/,
     3   ykow(10)/2.35e-2/, ykow(11)/2.46e-2/
      data ukow(1)/4.0/, ukow(2)/2.0/, ukow(3)/1.0/,
     1   ukow(4)/5.0e-1/, ukow(5)/2.5e-1/, ukow(6)/1.67e-1/,
     2   ukow(7)/1.25e-1/, ukow(8)/1.0e-1/, ukow(9)/8.33e-2/,
     3   ukow(10)/7.14e-2/, ukow(11)/6.25e-2/
      data yosb1(1)/8.44e-1/, yosb1(2)/9.08e-1/, yosb1(3)/9.32e-1/,
     1   yosb1(4)/9.36e-1/, yosb1(5)/9.25e-1/, yosb1(6)/9.08e-1/,
     2   yosb1(7)/8.81e-1/, yosb1(8)/8.50e-1/, yosb1(9)/8.18e-1/,
     3   yosb1(10)/7.84e-1/, yosb1(11)/7.51e-1/, yosb1(12)/7.18e-1/,
     4   yosb1(13)/6.85e-1/, yosb1(14)/6.58e-1/, yosb1(15)/6.28e-1/,
     5   yosb1(16)/6.03e-1/, yosb1(17)/5.80e-1/, yosb1(18)/5.58e-1/,
     6   yosb1(19)/5.38e-1/, yosb1(20)/5.22e-1/, yosb1(21)/5.06e-1/,
     7   yosb1(22)/4.90e-1/, yosb1(23)/4.78e-1/, yosb1(24)/4.67e-1/,
     8   yosb1(25)/4.57e-1/, yosb1(26)/4.48e-1/, yosb1(27)/4.38e-1/,
     9   yosb1(28)/4.31e-1/, yosb1(29)/4.24e-1/, yosb1(30)/4.20e-1/,
     a   yosb1(31)/4.14e-1/, yosb1(32)/4.11e-1/, yosb1(33)/4.06e-1/
      data yosb2(1)/1.366/, yosb2(2)/1.191/, yosb2(3)/1.112/,
     1   yosb2(4)/1.013/, yosb2(5)/9.91e-1/, yosb2(6)/8.85e-1/,
     2   yosb2(7)/8.31e-1/, yosb2(8)/8.47e-1/, yosb2(9)/7.86e-1/,
     3   yosb2(10)/7.25e-1/, yosb2(11)/7.46e-1/, yosb2(12)/6.79e-1/,
     4   yosb2(13)/6.08e-1/, yosb2(14)/6.55e-1/, yosb2(15)/6.16e-1/,
     5   yosb2(16)/6.06e-1/, yosb2(17)/6.02e-1/, yosb2(18)/6.26e-1/,
     6   yosb2(19)/6.51e-1/, yosb2(20)/7.24e-1/, yosb2(21)/6.49e-1/,
     7   yosb2(22)/6.49e-1/, yosb2(23)/6.94e-1/, yosb2(24)/6.44e-1/,
     8   yosb2(25)/6.24e-1/, yosb2(26)/6.61e-1/, yosb2(27)/6.12e-1/,
     9   yosb2(28)/5.58e-1/, yosb2(29)/5.33e-1/, yosb2(30)/4.95e-1/,
     a   yosb2(31)/5.00e-1/, yosb2(32)/4.23e-1/, yosb2(33)/3.95e-1/,
     b   yosb2(34)/3.75e-1/, yosb2(35)/3.72e-1/, yosb2(36)/3.91e-1/,
     c   yosb2(37)/3.96e-1/, yosb2(38)/4.05e-1/, yosb2(39)/4.28e-1/,
     d   yosb2(40)/4.29e-1/, yosb2(41)/5.23e-1/, yosb2(42)/5.62e-1/,
     e   yosb2(43)/6.07e-1/, yosb2(44)/6.53e-1/, yosb2(45)/6.72e-1/,
     f   yosb2(46)/7.08e-1/, yosb2(47)/6.33e-1/, yosb2(48)/6.68e-1/,
     g   yosb2(49)/6.45e-1/, yosb2(50)/6.32e-1/, yosb2(51)/5.91e-1/,
     h   yosb2(52)/5.59e-1/, yosb2(53)/5.97e-1/, yosb2(54)/6.25e-1/,
     i   yosb2(55)/7.39e-1/, yosb2(56)/7.10e-1/, yosb2(57)/7.29e-1/,
     j   yosb2(58)/7.20e-1/, yosb2(59)/6.36e-1/, yosb2(60)/5.81e-1/
      data yosb2(61)/4.28e-1/, yosb2(62)/2.92e-1/, yosb2(63)/1.62e-1/,
     1   yosb2(64)/9.8e-2/, yosb2(65)/5.4e-2/
      data ymeyer(1)/3.478e4/, ymeyer(2)/2.861e4/, ymeyer(3)/2.365e4/,
     1   ymeyer(4)/1.963e4/, ymeyer(5)/1.637e4/, ymeyer(6)/1.372e4/,
     2   ymeyer(7)/1.154e4/, ymeyer(8)/9.744e3/, ymeyer(9)/8.261e3/,
     3   ymeyer(10)/7.030e3/, ymeyer(11)/6.005e3/, ymeyer(12)/5.147e3/,
     4   ymeyer(13)/4.427e3/, ymeyer(14)/3.820e3/, ymeyer(15)/3.307e3/,
     5   ymeyer(16)/2.872e3/
c
      data expmax/0./, uftolg/0./
c
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
 100  r(1) = 1.0e1*(x(2) - x(1)**2)
      r(2) = 1.0 - x(1)
      go to 9999
c  ***  helix   ***
 200  theta = atan2(x(2), x(1))/twopi
      if (x(1) .le. 0. .and. x(2) .le. 0.) theta = theta + 1.
      r(1) = 1.0e1*(x(3) - 1.0e1*theta)
      r(2) = 1.0e1*(sqrt(x(1)**2 + x(2)**2) - 1.0)
      r(3) = x(3)
      go to 9999
c  ***  singular   ***
 300  r(1) = x(1) + 1.0e1*x(2)
      r(2) = sqrt(5.0)*(x(3) - x(4))
      r(3) = (x(2) - 2.0*x(3))**2
      r(4) = sqrt(1.0e1)*(x(1) - x(4))**2
      go to 9999
c  ***  woods   ***
 400  r(1) = 1.0e1*(x(2) - x(1)**2)
      r(2) = 1.0 - x(1)
      r(3) = sqrt(9.0e1)*(x(4) - x(3)**2)
      r(4) = 1.0 - x(3)
      r(5) = sqrt(9.9)*(x(2) + x(4) - 2.)
      t = sqrt(2.0e-1)
      r(6) = t*(x(2) - 1.0)
      r(7) = t*(x(4) - 1.0)
      go to 9999
c  ***  zangwill
 500  r(1) = x(1) - x(2) + x(3)
      r(2) = -x(1) + x(2) + x(3)
      r(3) = x(1) + x(2) - x(3)
      go to 9999
c  ***  engvall   ***
 600  r(1) = x(1)**2 + x(2)**2 + x(3)**2 - 1.0
      r(2) = x(1)**2 + x(2)**2 + (x(3) - 2.0)**2 - 1.0
      r(3) = x(1) + x(2) + x(3) - 1.0
      r(4) = x(1) + x(2) - x(3) + 1.0
      r(5) = x(1)**3 + 3.0*x(2)**2 + (5.0*x(3) - x(1) + 1.0)**2
     1               - 3.6e1
      go to 9999
c  ***  branin ***
 700  r(1) = 4.0*(x(1) + x(2))
      r(2) = r(1) + (x(1) - x(2))*((x(1) - 2.0)**2 +
     1       x(2)**2 - 1.0)
      go to 9999
c  ***  beale  ***
 800  r(1) = 1.5 - x(1)*(1.0 - x(2))
      r(2) = 2.25 - x(1)*(1.0 - x(2)**2)
      r(3) = 2.625 - x(1)*(1.0 -  x(2)**3)
      go to 9999
c  ***  cragg and levy  ***
 900  r(1) = (exp(x(1)) - x(2))**2
      r(2) = 1.0e1*(x(2) - x(3))**3
      r(3) = ( sin(x(3) - x(4)) / cos(x(3) - x(4)) )**2
      r(4) = x(1)**4
      r(5) = x(4) - 1.0
      go to 9999
c  ***  box  ***
 1000 if (expmax .gt. 0.) go to 1001
         expmax = 1.999 * alog(rmdcon(5))
         expmin = 1.999 * alog(rmdcon(2))
 1001 if (-expmax .ge. amin1(x(1), x(2), x(3))) go to 1003
      do 1002 i = 1,10
         ti = -0.1*float(i)
         t1 = ti*x(1)
         e1 = 0.
         if (t1 .gt. expmin) e1 = exp(t1)
         t2 = ti*x(2)
         e2 = 0.
         if (t2 .gt. expmin) e2 = exp(t2)
         r(i) = (e1 - e2) - x(3)*(exp(ti) - exp(1.0e1*ti))
 1002 continue
      go to 9999
 1003 nfcall = -1
      go to 9999
c  ***  davidon 1  ***
 1100 nm1 = n - 1
      do 1102 i = 1, nm1
         r1 = 0.0
         ti = float(i)
         t = 1.
         do 1101 j = 1,p
              r1 = r1 + t*x(j)
              t = t*ti
 1101         continue
         r(i) = r1
 1102    continue
      r(n) = x(1) - 1.0
      go to 9999
c  ***  freudenstein and roth  ***
 1200 r(1) = -1.3e1 + x(1) - 2.0*x(2) + 5.0*x(2)**2 - x(2)**3
      r(2) = -2.9e1 + x(1) - 1.4e1*x(2) + x(2)**2 + x(2)**3
      go to 9999
c  ***  watson  ***
 1300  continue
 1400  continue
 1500  continue
 1600 do 1602 i = 1, 29
         ti = float(i)/2.9e1
         r1 = 0.0
         r2 = x(1)
         t = 1.0
         do 1601 j = 2, p
              r1 = r1 + float(j-1)*t*x(j)
              t = t*ti
              r2 = r2 + t*x(j)
 1601         continue
         r(i) = r1 - r2*r2 - 1.0
         if (nex .ge. 33 .and. nex .le. 36) r(i) = r(i) + 10.
 1602    continue
      r(30) = x(1)
      r(31) = x(2) - x(1)**2 - 1.0
      if (nex .lt. 33 .or. nex .gt. 36) go to 9999
      r(30) = r(30) + 10.
      r(31) = r(31) + 10.
      go to 9999
c  ***  chebyquad  ***
 1700 do 1701 i = 1,n
 1701    r(i) = 0.0
      do 1702 j = 1,n
         tim1 = 1.0
         ti = 2.0*x(j) - 1.0
         z = ti + ti
         do 1702 i = 1,n
              r(i) = r(i) + ti
              tip1 = z*ti -tim1
              tim1 = ti
              ti = tip1
 1702         continue
      floatn = float(n)
      do 1703 i = 1,n
         ti = 0.0
         if (mod(i,2) .eq. 0) ti = -1.0/float(i*i - 1)
         r(i) = ti - r(i)/floatn
 1703    continue
      go to 9999
c  ***  brown and dennis  ***
 1800  do 1801 i = 1, n
         ti = 0.2*float(i)
         r(i) = (x(1) + x(2)*ti - exp(ti))**2 +
     1             (x(3) + x(4)*sin(ti) - cos(ti))**2
 1801    continue
      go to 9999
c  ***  bard  ***
 1900 do 1901 i = 1, 15
         u = float(i)
         v = 1.6e1 - u
         w = amin1(u,v)
         r(i) = ybard(i) - (x(1) + u/(x(2)*v + x(3)*w))
         if (nex .eq. 30) r(i) = r(i) + 10.
 1901    continue
      go to 9999
c  ***  jennrich and sampson  ***
 2000 do 2001 i = 1, 10
         ti = float(i)
         r(i) = 2.0 + 2.0*ti - (exp(ti*x(1)) +
     1          exp(ti*x(2)))
 2001    continue
      go to 9999
c  ***  kowalik and osborne  ***
 2100 do 2101 i = 1, 11
         r(i) = ykow(i) - x(1)*(ukow(i)**2 + x(2)*ukow(i))/(ukow(i)**2 +
     1          x(3)*ukow(i) + x(4))
         if (nex .eq. 31) r(i) = r(i) + 10.
 2101    continue
      go to 9999
c  ***  osborne 1  ***
 2200 do 2201 i = 1, 33
         ti = 1.0e1*float(1-i)
         r(i) = yosb1(i) - (x(1) + x(2)*exp(x(4)*ti) +
     1          x(3)*exp(x(5)*ti))
 2201    continue
      go to 9999
c  ***  osborne 2  ***
c     ***  uftolg is a machine-dependent constant.  it is just slightly
c     ***  larger than the log of the smallest positive machine number.
 2300 if (uftolg .eq. 0.) uftolg = 1.999 * alog(rmdcon(2))
      do 2302 i = 1, 65
         ti = 0.1*float(1-i)
         ri = x(1)*exp(x(5)*ti)
         do 2301 j = 2, 4
              t = 0.
              theta = -x(j+4) * (ti + x(j+7))**2
              if (theta .gt. uftolg) t = exp(theta)
              ri = ri + x(j)*t
 2301         continue
         r(i) = yosb2(i) - ri
 2302 continue
      go to 9999
c  ***  madsen  ***
 2400 r(1) = x(1)**2 + x(2)**2 + x(1)*x(2)
      r(2) = sin(x(1))
      r(3) = cos(x(2))
      go to 9999
c  ***  meyer  ***
 2500 do 2501 i = 1, 16
         ti = float(5*i + 45)
         r(i)=x(1)*exp(x(2)/(ti + x(3))) - ymeyer(i)
         if (nex .eq. 32) r(i) = r(i) + 10.
 2501    continue
      go to 9999
c  ***  brown  ***
 2600 continue
 2700 continue
 2800 continue
 2900 t = x(1) - float(n + 1)
      do 2901 i = 2, n
 2901    t = t + x(i)
      nm1 = n - 1
      do 2902 i = 1, nm1
 2902    r(i) = t + x(i)
      t = x(1)
      do 2903 i = 2, n
 2903    t = t * x(i)
      r(n) = t - 1.0
      go to 9999
c
 9999 return
c     ..... last card of testr .........................................
      end
      subroutine today(datime)
c
c  ***  supply sumsol version  ***
c
c/6
      real datime(4), dt1, dt2, dt3, dt4
      data dt1,dt2,dt3,dt4/4hnl2s,4hol  ,4hver.,4h2.2 /
c/7
c     character*4 datime(4), dt1, dt2, dt3, dt4
c     data dt1,dt2,dt3,dt4/'nl2s','ol  ','ver.','2.2 '/
c/
c
      datime(1) = dt1
      datime(2) = dt2
      datime(3) = dt3
      datime(4) = dt4
 999  return
c  ***  last line of datime follows  ***
      end
      subroutine xinit(p, x, nex)
c
c     *****parameters...
c
      integer nex, p
      real x(p)
c
c     ..................................................................
c
c     *****purpose...
c     this routine initializes the solution vector x according to
c     the initial values for the various test functions given in
c     references (1), (2), and (3).
c     subroutines testr and testj.  (see testr for references.)
c
c     *****parameter description...
c     on input...
c
c        nex is the test problem number.
c
c        p is the number of parameters.
c
c     on output...
c
c        x is the initial guess to the solution.
c
c     *****application and usage restrictions...
c     this routine is called by nltest.
c
c     ..................................................................
c
c     *****local variables...
      integer i
      real pp1inv
c     *****intrinsic functions...
c/+
      real float
c/
c
      go to (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100,
     1   1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100,
     2   2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 1900, 2100,
     3   2500, 1300, 1400, 1500, 1600),nex
c
c  ***  rosenbrock  ***
 100  x(1) = -1.2
      x(2) = 1.0
      go to 9999
c  ***  helix  ***
 200  x(1) = -1.0
      x(2) = 0.0
      x(3) = 0.0
      go to 9999
c  *** singular  ***
 300  x(1) = 3.0
      x(2) = -1.0
      x(3) = 0.0
      x(4) = 1.0
      go to 9999
c  ***  woods  ***
 400  x(1) = -3.0
      x(2) = -1.0
      x(3) = -3.0
      x(4) = -1.0
      go to 9999
c  ***  zangwill  ***
 500  x(1) = 1.0e2
      x(2) = -1.0
      x(3) = 2.5
      go to 9999
c  ***  engvall  ***
 600  x(1) = 1.0
      x(2) = 2.0
      x(3) = 0.0
      go to 9999
c  *** branin  ***
 700  x(1) = 2.0
      x(2) = 0.0
      go to 9999
c  ***  beale  ***
 800  x(1) = 1.0e-1
      x(2) = 1.0e-1
      go to 9999
c  *** cragg and levy  ***
 900  x(1) = 1.0
      x(2) = 2.0
      x(3) = 2.0
      x(4) = 2.0
      go to 9999
c  ***  box  ***
 1000 x(1) = 0.0
      x(2) = 1.0e1
      x(3) = 2.0e1
      go to 9999
c  ***  davidon 1  ***
 1100 do 1101 i = 1,p
 1101    x(i) = 0.0
      go to 9999
c  ***  freudenstein and roth  ***
 1200 x(1) = 1.5e1
      x(2) = -2.0
      go to 9999
c  ***  watson  ***
 1300 continue
 1400 continue
 1500 continue
 1600 do 1601 i = 1,p
 1601    x(i) = 0.0
      go to 9999
c  ***  chebyquad  ***
 1700 pp1inv = 1.0/float(p + 1)
      do 1701 i = 1, p
 1701    x(i) = float(i)*pp1inv
      go to 9999
c  *** brown and dennis  ***
 1800 x(1) = 2.5e1
      x(2) = 5.0
      x(3) = -5.0
      x(4) = -1.0
      go to 9999
c  ***  bard  ***
 1900 x(1) = 1.
      x(2) = 1.
      x(3) = 1.
      go to 9999
c  ***  jennrich and sampson  ***
 2000 x(1) = 3.0e-1
      x(2) = 4.0e-1
      go to 9999
c  ***  kowalik and osborne  ***
 2100 x(1) = 2.5e-1
      x(2) = 3.9e-1
      x(3) = 4.15e-1
      x(4) = 3.9e-1
      go to 9999
c  ***  osborne 1  ***
 2200 x(1) = 5.0e-1
      x(2) = 1.5
      x(3) = -1.0
      x(4) = 1.0e-2
      x(5) = 2.0e-2
      go to 9999
c  ***  osborne 2  ***
 2300 x(1) = 1.3
      x(2) = 6.5e-1
      x(3) = 6.5e-1
      x(4) = 7.0e-1
      x(5) = 6.0e-1
      x(6) = 3.0
      x(7) = 5.0
      x(8) = 7.0
      x(9) = 2.0
      x(10) = 4.5
      x(11) = 5.5
      go to 9999
c  ***  madsen  ***
 2400 x(1) = 3.0
      x(2) = 1.0
      go to 9999
c  ***  meyer  **
 2500 x(1) = 2.0e-2
      x(2) = 4.0e3
      x(3) = 2.5e2
      go to 9999
c  ***  brown  ***
 2600 continue
 2700 continue
 2800 continue
 2900 do 2901 i = 1, p
 2901    x(i) = 5.e-1
      go to 9999
c
c
 9999 return
      end
c///////////////////////////////////////////////////////////////////////
c     ***  test nl2sol and nl2sno on madsen example  ***
      integer iv(62), uiparm(1)
      double precision v(147), x(2), urparm(1)
      external madr, madj
      x(1) = 3.0d0
      x(2) = 1.0d0
      iv(1) = 0
      call nl2sol(3, 2, x, madr, madj, iv, v, uiparm, urparm, madr)
      iv(1) = 12
      x(1) = 3.0d0
      x(2) = 1.0d0
      call nl2sno(3, 2, x, madr, iv, v, uiparm, urparm, madr)
      stop
      end
      subroutine madr(n, p, x, nf, r, uiparm, urparm, ufparm)
      integer n, p, nf, uiparm(1)
      double precision x(p), r(n), urparm(1)
      external ufparm
      r(1) = x(1)**2 + x(2)**2 + x(1)*x(2)
      r(2) = dsin(x(1))
      r(3) = dcos(x(2))
      return
      end
      subroutine madj(n, p, x, nf, j, uiparm, urparm, ufparm)
      integer n, p, nf, uiparm(1)
      double precision x(p), j(n,p), urparm(1)
      external ufparm
      j(1,1) = 2.0d0*x(1) + x(2)
      j(1,2) = 2.0d0*x(2) + x(1)
      j(2,1) = dcos(x(1))
      j(2,2) = 0.0d0
      j(3,1) = 0.0d0
      j(3,2) = -dsin(x(2))
      return
      end
c///////////////////////////////////////////////////////////////////////
      integer function imdcon(k)
c
      integer k
c
c  ***  return integer machine-dependent constants  ***
c
c     ***  k = 1 means return standard output unit number.   ***
c     ***  k = 2 means return alternate output unit number.  ***
c     ***  k = 3 means return  input unit number.            ***
c          (note -- k = 2, 3 are used only by test programs.)
c
      integer mdcon(3)
      data mdcon(1)/6/, mdcon(2)/8/, mdcon(3)/5/
c
      imdcon = mdcon(k)
      return
c  ***  last card of imdcon follows  ***
      end
      double precision function rmdcon(k)
c
c  ***  return machine dependent constants used by nl2sol  ***
c
c +++  comments below contain data statements for various machines.  +++
c +++  to convert to another machine, place a c in column 1 of the   +++
c +++  data statement line(s) that correspond to the current machine +++
c +++  and remove the c from column 1 of the data statement line(s)  +++
c +++  that correspond to the new machine.                           +++
c
      integer k
c
c  ***  the constant returned depends on k...
c
c  ***        k = 1... smallest pos. eta such that -eta exists.
c  ***        k = 2... square root of 1.001*eta.
c  ***        k = 3... unit roundoff = smallest pos. no. machep such
c  ***                 that 1 + machep .gt. 1 .and. 1 - machep .lt. 1.
c  ***        k = 4... square root of 0.999*machep.
c  ***        k = 5... square root of 0.999*big (see k = 6).
c  ***        k = 6... largest machine no. big such that -big exists.
c
      double precision big, eta, machep
c/+
      double precision dsqrt
c/
      double precision one001, pt999
c
      data one001/1.001d0/, pt999/0.999d0/
c
c  +++  ibm 360, ibm 370, or xerox  +++
c
      data big/z7fffffffffffffff/, eta/z0010000000000000/,
     1     machep/z3410000000000000/
c
c  +++  data general  +++
c
c     data big/0.7237005577d+76/, eta/0.5397605347d-78/,
c    1     machep/2.22044605d-16/
c
c  +++  dec 11  +++
c
c     data big/1.7d+38/, eta/2.938735878d-39/, machep/2.775557562d-17/
c
c  +++  hp3000  +++
c
c     data big/1.157920892d+77/, eta/8.636168556d-78/,
c    1     machep/5.551115124d-17/
c
c  +++  honeywell  +++
c
c     data big/1.69d+38/, eta/5.9d-39/, machep/2.1680435d-19/
c
c  +++  dec10  +++
c
c     data big/"377777100000000000000000/,
c    1     eta/"002400400000000000000000/,
c    2     machep/"104400000000000000000000/
c
c  +++  burroughs  +++
c
c     data big/o0777777777777777,o7777777777777777/,
c    1     eta/o1771000000000000,o7770000000000000/,
c    2     machep/o1451000000000000,o0000000000000000/
c
c  +++  control data  +++
c
c
c     data big/37767777777777777777b,37167777777777777777b/,
c    1     eta/00014000000000000000b,00000000000000000000b/,
c    2     machep/15614000000000000000b,15010000000000000000b/
c
c  +++  prime  +++
c
c     data big/1.0d+9786/, eta/1.0d-9860/, machep/1.4210855d-14/
c
c  +++  univac  +++
c
c     data big/8.988d+307/, eta/1.2d-308/, machep/1.734723476d-18/
c
c  +++  vax  +++
c
c     data big/1.7d+38/, eta/2.939d-39/, machep/1.3877788d-17/
c
c-------------------------------  body  --------------------------------
c
      go to (10, 20, 30, 40, 50, 60), k
c
 10   rmdcon = eta
      go to 999
c
 20   rmdcon = dsqrt(one001*eta)
      go to 999
c
 30   rmdcon = machep
      go to 999
c
 40   rmdcon = dsqrt(pt999*machep)
      go to 999
c
 50   rmdcon = dsqrt(pt999*big)
      go to 999
c
 60   rmdcon = big
c
 999  return
c  ***  last card of rmdcon follows  ***
      end
c///////////////////////////////////////////////////////////////////////
      subroutine nl2sol(n, p, x, calcr, calcj, iv, v, uiparm, urparm,
     1                  ufparm)
c
c  ***  minimize nonlinear sum of squares using analytic jacobian  ***
c  ***  (nl2sol version 2.2)  ***
c
      integer n, p, iv(1), uiparm(1)
      double precision x(p), v(1), urparm(1)
c     dimension iv(60+p),  v(93 + n*p + 3*n + p*(3*p+33)/2)
c     dimension uiparm(*), urparm(*)
      external calcr, calcj, ufparm
c
c  ***  purpose  ***
c
c        given a p-vector x of parameters, calcr computes an n-vector
c     r = r(x) of residuals corresponding to x.  (r(x) probably arises
c     from a nonlinear model involving p parameters and n observations.)
c     this routine interacts with nl2itr to seek a parameter vector x
c     that minimizes the sum of the squares of (the components of) r(x),
c     i.e., that minimizes the sum-of-squares function
c     f(x) = (r(x)**t) * r(x) / 2.  r(x) is assumed to be a twice con-
c     tinuously differentiable function of x.
c
c--------------------------  parameter usage  --------------------------
c
c n........ (input) the number of observations, i.e., the number of
c                  components in r(x).  n must be .ge. p.
c p........ (input) the number of parameters (components in x).  p must
c                  be positive.
c x........ (input/output).  on input, x is an initial guess at the
c                  desired parameter estimate.  on output, x contains
c                  the best parameter estimate found.
c calcr.... (input) a subroutine which, given x, computes r(x).  calcr
c                  must be declared external in the calling program.
c                  it is invoked by
c                       call calcr(n,p,x,nf,r,uiparm,urparm,ufparm)
c                  when calcr is called, nf is the invocation count
c                  for calcr.  it is included for possible use with
c                  calcj.  if x is out of bounds (e.g. if it would
c                  cause overflow in computing r(x)), then calcr should
c                  set nf to 0.  this will cause a shorter step to be
c                  attempted.  the other parameters are as described
c                  above and below.  calcr should not change n, p, or x.
c calcj.... (input) a subroutine which, given x, computes the jacobian
c                  matrix j of r at x, i.e., the n by p matrix whose
c                  (i,k) entry is the partial derivative of the i-th
c                  component of r with respect to x(k).  calcj must be
c                  declared external in the calling program.  it is
c                  invoked by
c                       call calcj(n,p,x,nf,j,uiparm,urparm,ufparm)
c                  nf is the invocation count for calcr at the time
c                  r(x) was evaluated.  the x passed to calcj is
c                  usually the one passed to calcr on either its most
c                  recent invocation or the one prior to it.  if calcr
c                  saves intermediate results for use by calcj, then it
c                  is possible to tell from nf whether they are valid
c                  for the current x (or which copy is valid if two
c                  copies are kept).  if j cannot be computed at x,
c                  then calcj should set nf to 0.  in this case, nl2sol
c                  will return with iv(1) = 15.  the other parameters
c                  to calcj are as described above and below.  calcj
c                  should not change n, p, or x.
c iv....... (input/output) an integer value array of length at least
c                  60 + p that helps control the nl2sol algorithm and
c                  that is used to store various intermediate quanti-
c                  ties.  of particular interest are the initialization/
c                  return code iv(1) and the entries in iv that control
c                  printing and limit the number of iterations and func-
c                  tion evaluations.  see the section on iv input
c                  values below.
c v........ (input/output) a floating-point value array of length at
c                  least 93 + n*p + 3*n + p*(3*p+33)/2 that helps con-
c                  trol the nl2sol algorithm and that is used to store
c                  various intermediate quantities.  of particular in-
c                  terest are the entries in v that limit the length of
c                  the first step attempted (lmax0), specify conver-
c                  gence tolerances (afctol, rfctol, xctol, xftol),
c                  and help choose the step size used in computing the
c                  covariance matrix (delta0).  see the section on
c                  (selected) v input values below.
c uiparm... (input) user integer parameter array passed without change
c                  to calcr and calcj.
c urparm... (input) user floating-point parameter array passed without
c                  change to calcr and calcj.
c ufparm... (input) user external subroutine or function passed without
c                  change to calcr and calcj.
c
c  ***  iv input values (from subroutine dfault)  ***
c
c iv(1)...  on input, iv(1) should have a value between 0 and 12......
c             0 and 12 mean this is a fresh start.  0 means that
c             dfault(iv, v) is to be called to provide all default
c             values to iv and v.  12 (the value that dfault assigns to
c             iv(1)) means the caller has already called dfault(iv, v)
c             and has possibly changed some iv and/or v entries to non-
c             default values.  default = 12.
c iv(covprt)... iv(14) = 1 means print a covariance matrix at the solu-
c             tion.  (this matrix is computed just before a return with
c             iv(1) = 3, 4, 5, 6.)
c             iv(covprt) = 0 means skip this printing.  default = 1.
c iv(covreq)... iv(15) = nonzero means compute a covariance matrix
c             just before a return with iv(1) = 3, 4, 5, 6.  in
c             this case, an approximate covariance matrix is obtained
c             in one of several ways.  let k = abs(iv(covreq)) and let
c             scale = 2*f(x)/max(1,n-p),  where 2*f(x) is the residual
c             sum of squares.  if k = 1 or 2, then a finite-difference
c             hessian approximation h is obtained.  if h is positive
c             definite (or, for k = 3, if the jacobian matrix j at x
c             is nonsingular), then one of the following is computed...
c                  k = 1....  scale * h**-1 * (j**t * j) * h**-1.
c                  k = 2....  scale * h**-1.
c                  k = 3....  scale * (j**t * j)**-1.
c             (j**t is the transpose of j, while **-1 means inverse.)
c             if iv(covreq) is positive, then both function and grad-
c             ient values (calls on calcr and calcj) are used in com-
c             puting h (with step sizes determined using v(delta0) --
c             see below), while if iv(covreq) is negative, then only
c             function values (calls on calcr) are used (with step
c             sizes determined using v(dltfdc) -- see below).  if
c             iv(covreq) = 0, then no attempt is made to compute a co-
c             variance matrix (unless iv(covprt) = 1, in which case
c             iv(covreq) = 1 is assumed).  see iv(covmat) below.
c             default = 1.
c iv(dtype).... iv(16) tells how the scale vector d (see ref. 1) should
c             be chosen.  iv(dtype) .ge. 1 means choose d as described
c             below with v(dfac).  iv(dtype) .le. 0 means the caller
c             has chosen d and has stored it in v starting at
c             v(94 + 2*n + p*(3*p + 31)/2).  default = 1.
c iv(inits).... iv(25) tells how the s matrix (see ref. 1) should be
c             initialized.  0 means initialize s to 0 (and start with
c             the gauss-newton model).  1 and 2 mean that the caller
c             has stored the lower triangle of the initial s rowwise in
c             v starting at v(87+2*p).  iv(inits) = 1 means start with
c             the gauss-newton model, while iv(inits) = 2 means start
c             with the augmented model (see ref. 1).  default = 0.
c iv(mxfcal)... iv(17) gives the maximum number of function evaluations
c             (calls on calcr, excluding those used to compute the co-
c             variance matrix) allowed.  if this number does not suf-
c             fice, then nl2sol returns with iv(1) = 9.  default = 200.
c iv(mxiter)... iv(18) gives the maximum number of iterations allowed.
c             it also indirectly limits the number of gradient evalua-
c             tions (calls on calcj, excluding those used to compute
c             the covariance matrix) to iv(mxiter) + 1.  if iv(mxiter)
c             iterations do not suffice, then nl2sol returns with
c             iv(1) = 10.  default = 150.
c iv(outlev)... iv(19) controls the number and length of iteration sum-
c             mary lines printed (by itsmry).  iv(outlev) = 0 means do
c             not print any summary lines.  otherwise, print a summary
c             line after each abs(iv(outlev)) iterations.  if iv(outlev)
c             is positive, then summary lines of length 117 (plus carri-
c             age control) are printed, including the following...  the
c             iteration and function evaluation counts, current func-
c             tion value (v(f) = half the sum of squares), relative
c             difference in function values achieved by the latest step
c             (i.e., reldf = (f0-v(f))/f0, where f0 is the function
c             value from the previous iteration), the relative function
c             reduction predicted for the step just taken (i.e.,
c             preldf = v(preduc) / f0, where v(preduc) is described
c             below), the scaled relative change in x (see v(reldx)
c             below), the models used in the current iteration (g =
c             gauss-newton, s=augmented), the marquardt parameter
c             stppar used in computing the last step, the sizing factor
c             used in updating s, the 2-norm of the scale vector d
c             times the step just taken (see ref. 1), and npreldf, i.e.,
c             v(nreduc)/f0, where v(nreduc) is described below -- if
c             npreldf is positive, then it is the relative function
c             reduction predicted for a newton step (one with
c             stppar = 0).  if npreldf is zero, either the gradient
c             vanishes (as does preldf) or else the augmented model
c             is being used and its hessian is indefinite (with preldf
c             positive).  if npreldf is negative, then it is the nega-
c             of the relative function reduction predicted for a step
c             computed with step bound v(lmax0) for use in testing for
c             singular convergence.
c                  if iv(outlev) is negative, then lines of maximum
c             length 79 (or 55 is iv(covprt) = 0) are printed, includ-
c             ing only the first 6 items listed above (through reldx).
c             default = 1.
c iv(parprt)... iv(20) = 1 means print any nondefault v values on a
c             fresh start or any changed v values on a restart.
c             iv(parprt) = 0 means skip this printing.  default = 1.
c iv(prunit)... iv(21) is the output unit number on which all printing
c             is done.  iv(prunit) = 0 means suppress all printing.
c             (setting iv(prunit) to 0 is the only way to suppress the
c             one-line termination reason message printed by itsmry.)
c             default = standard output unit (unit 6 on most systems).
c iv(solprt)... iv(22) = 1 means print out the value of x returned (as
c             well as the corresponding gradient and scale vector d).
c             iv(solprt) = 0 means skip this printing.  default = 1.
c iv(statpr)... iv(23) = 1 means print summary statistics upon return-
c             ing.  these consist of the function value (half the sum
c             of squares) at x, v(reldx) (see below), the number of
c             function and gradient evaluations (calls on calcr and
c             calcj respectively, excluding any calls used to compute
c             the covariance), the relative function reductions predict-
c             ed for the last step taken and for a newton step (or per-
c             haps a step bounded by v(lmax0) -- see the descriptions
c             of preldf and npreldf under iv(outlev) above), and (if an
c             attempt was made to compute the covariance) the number of
c             calls on calcr and calcj used in trying to compute the
c             covariance.  iv(statpr) = 0 means skip this printing.
c             default = 1.
c iv(x0prt).... iv(24) = 1 means print the initial x and scale vector d
c             (on a fresh start only).  iv(x0prt) = 0 means skip this
c             printing.  default = 1.
c
c  ***  (selected) iv output values  ***
c
c iv(1)........ on output, iv(1) is a return code....
c             3 = x-convergence.  the scaled relative difference be-
c                  tween the current parameter vector x and a locally
c                  optimal parameter vector is very likely at most
c                  v(xctol).
c             4 = relative function convergence.  the relative differ-
c                  ence between the current function value and its lo-
c                  cally optimal value is very likely at most v(rfctol).
c             5 = both x- and relative function convergence (i.e., the
c                  conditions for iv(1) = 3 and iv(1) = 4 both hold).
c             6 = absolute function convergence.  the current function
c                  value is at most v(afctol) in absolute value.
c             7 = singular convergence.  the hessian near the current
c                  iterate appears to be singular or nearly so, and a
c                  step of length at most v(lmax0) is unlikely to yield
c                  a relative function decrease of more than v(rfctol).
c             8 = false convergence.  the iterates appear to be converg-
c                  ing to a noncritical point.  this may mean that the
c                  convergence tolerances (v(afctol), v(rfctol),
c                  v(xctol)) are too small for the accuracy to which
c                  the function and gradient are being computed, that
c                  there is an error in computing the gradient, or that
c                  the function or gradient is discontinuous near x.
c             9 = function evaluation limit reached without other con-
c                  vergence (see iv(mxfcal)).
c            10 = iteration limit reached without other convergence
c                  (see iv(mxiter)).
c            11 = stopx returned .true. (external interrupt).  see the
c                  usage notes below.
c            13 = f(x) cannot be computed at the initial x.
c            14 = bad parameters passed to assess (which should not
c                  occur).
c            15 = the jacobian could not be computed at x (see calcj
c                  above).
c            16 = n or p (or parameter nn to nl2itr) out of range --
c                  p .le. 0 or n .lt. p or nn .lt. n.
c            17 = restart attempted with n or p (or par. nn to nl2itr)
c                  changed.
c            18 = iv(inits) is out of range.
c            19...45 = v(iv(1)) is out of range.
c            50 = iv(1) was out of range.
c            87...(86+p) = jtol(iv(1)-86) (i.e., v(iv(1)) is not
c                  positive (see v(dfac) below).
c iv(covmat)... iv(26) tells whether a covariance matrix was computed.
c             if (iv(covmat) is positive, then the lower triangle of
c             the covariance matrix is stored rowwise in v starting at
c             v(iv(covmat)).  if iv(covmat) = 0, then no attempt was
c             made to compute the covariance.  if iv(covmat) = -1,
c             then the finite-difference hessian was indefinite.  and
c             and if iv(covmat) = -2, then a successful finite-differ-
c             encing step could not be found for some component of x
c             (i.e., calcr set nf to 0 for each of two trial steps).
c             note that iv(covmat) is reset to 0 after each successful
c             step, so if such a step is taken after a restart, then
c             the covariance matrix will be recomputed.
c iv(d)........ iv(27) is the starting subscript in v of the current
c             scale vector d.
c iv(g)........ iv(28) is the starting subscript in v of the current
c             least-squares gradient vector (j**t)*r.
c iv(nfcall)... iv(6) is the number of calls so far made on calcr (i.e.,
c             function evaluations, including those used in computing
c             the covariance).
c iv(nfcov).... iv(40) is the number of calls made on calcr when
c             trying to compute covariance matrices.
c iv(ngcall)... iv(30) is the number of gradient evaluations (calls on
c             calcj) so far done (including those used for computing
c             the covariance).
c iv(ngcov).... iv(41) is the number of calls made on calcj when
c             trying to compute covariance matrices.
c iv(niter).... iv(31) is the number of iterations performed.
c iv(r)........ iv(50) is the starting subscript in v of the residual
c             vector r corresponding to x.
c
c  ***  (selected) v input values (from subroutine dfault)  ***
c
c v(afctol)... v(31) is the absolute function convergence tolerance.
c             if nl2sol finds a point where the function value (half
c             the sum of squares) is less than v(afctol), and if nl2sol
c             does not return with iv(1) = 3, 4, or 5, then it returns
c             with iv(1) = 6.  default = max(10**-20, machep**2), where
c             machep is the unit roundoff.
c v(delta0)... v(44) is a factor used in choosing the finite-difference
c             step size used in computing the covariance matrix when
c             iv(covreq) = 1 or 2.  for component i, step size
c                  v(delta0) * max(abs(x(i)), 1/d(i)) * sign(x(i))
c             is used, where d is the current scale vector (see ref. 1).
c             (if this step results in calcr setting nf to 0, then -0.5
c             times this step is also tried.)  default = machep**0.5,
c             where machep is the unit roundoff.
c v(dfac)..... v(41) and the d0 and jtol arrays (see v(d0init) and
c             v(jtinit)) are used in updating the scale vector d when
c             iv(dtype) .gt. 0.  (d is initialized according to
c             v(dinit).)  let d1(i) =
c               max(sqrt(jcnorm(i)**2 + max(s(i,i),0)), v(dfac)*d(i)),
c             where jcnorm(i) is the 2-norm of the i-th column of the
c             current jacobian matrix and s is the s matrix of ref. 1.
c             if iv(dtype) = 1, then d(i) is set to d1(i) unless
c             d1(i) .lt. jtol(i), in which case d(i) is set to
c                                max(d0(i), jtol(i)).
c             if iv(dtype) .ge. 2, then d is updated during the first
c             iteration as for iv(dtype) = 1 (after any initialization
c             due to v(dinit)) and is left unchanged thereafter.
c             default = 0.6.
c v(dinit).... v(38), if nonnegative, is the value to which the scale
c             vector d is initialized.  default = 0.
c v(dltfdc)... v(40) helps choose the step size used when computing the
c             covariance matrix when iv(covreq) = -1 or -2.  for
c             differences involving x(i), the step size first tried is
c                       v(dltfdc) * max(abs(x(i)), 1/d(i)),
c             where d is the current scale vector (see ref. 1).  (if
c             this step is too big the first time it is tried, i.e., if
c             calcr sets nf to 0, then -0.5 times this step is also
c             tried.)  default = machep**(1/3), where machep is the
c             unit roundoff.
c v(d0init)... v(37), if positive, is the value to which all components
c             of the d0 vector (see v(dfac)) are initialized.  if
c             v(dfac) = 0, then it is assumed that the caller has
c             stored d0 in v starting at v(p+87).  default = 1.0.
c v(jtinit)... v(39), if positive, is the value to which all components
c             of the jtol array (see v(dfac)) are initialized.  if
c             v(jtinit) = 0, then it is assumed that the caller has
c             stored jtol in v starting at v(87).  default = 10**-6.
c v(lmax0).... v(35) gives the maximum 2-norm allowed for d times the
c             very first step that nl2sol attempts.  it is also used
c             in testing for singular convergence -- if the function
c             reduction predicted for a step of length bounded by
c             v(lmax0) is at most v(rfctol) * abs(f0), where  f0  is
c             the function value at the start of the current iteration,
c             and if nl2sol does not return with iv(1) = 3, 4, 5, or 6,
c             then it returns with iv(1) = 7.    default = 100.
c v(rfctol)... v(32) is the relative function convergence tolerance.
c             if the current model predicts a maximum possible function
c             reduction (see v(nreduc)) of at most v(rfctol)*abs(f0) at
c             the start of the current iteration, where  f0  is the
c             then current function value, and if the last step attempt-
c             ed achieved no more than twice the predicted function
c             decrease, then nl2sol returns with iv(1) = 4 (or 5).
c             default = max(10**-10, machep**(2/3)), where machep is
c             the unit roundoff.
c v(tuner1)... v(26) helps decide when to check for false convergence
c             and to consider switching models.  this is done if the
c             actual function decrease from the current step is no more
c             than v(tuner1) times its predicted value.  default = 0.1.
c v(xctol).... v(33) is the x-convergence tolerance.  if a newton step
c             (see v(nreduc)) is tried that has v(reldx) .le. v(xctol)
c             and if this step yields at most twice the predicted func-
c             tion decrease, then nl2sol returns with iv(1) = 3 (or 5).
c             (see the description of v(reldx) below.)
c             default = machep**0.5, where machep is the unit roundoff.
c v(xftol).... v(34) is the false convergence tolerance.  if a step is
c             tried that gives no more than v(tuner1) times the predict-
c             ed function decrease and that has v(reldx) .le. v(xftol),
c             and if nl2sol does not return with iv(1) = 3, 4, 5, 6, or
c             7, then it returns with iv(1) = 8.  (see the description
c             of v(reldx) below.)  default = 100*machep, where
c             machep is the unit roundoff.
c v(*)........ dfault supplies to v a number of tuning constants, with
c             which it should ordinarily be unnecessary to tinker.  see
c             version 2.2 of the nl2sol usage summary (which is an
c             appendix to ref. 1).
c
c  ***  (selected) v output values  ***
c
c v(dgnorm)... v(1) is the 2-norm of (d**-1)*g, where g is the most re-
c             cently computed gradient and d is the corresponding scale
c             vector.
c v(dstnrm)... v(2) is the 2-norm of d*step, where step is the most re-
c             cently computed step and d is the current scale vector.
c v(f)........ v(10) is the current function value (half the sum of
c             squares).
c v(f0)....... v(13) is the function value at the start of the current
c             iteration.
c v(nreduc)... v(6), if positive, is the maximum function reduction
c             possible according to the current model, i.e., the func-
c             tion reduction predicted for a newton step (i.e.,
c             step = -h**-1 * g,  where  g = (j**t) * r  is the current
c             gradient and h is the current hessian approximation --
c             h = (j**t)*j  for the gauss-newton model and
c             h = (j**t)*j + s  for the augmented model).
c                  v(nreduc) = zero means h is not positive definite.
c                  if v(nreduc) is negative, then it is the negative of
c             the function reduction predicted for a step computed with
c             a step bound of v(lmax0) for use in testing for singular
c             convergence.
c v(preduc)... v(7) is the function reduction predicted (by the current
c             quadratic model) for the current step.  this (divided by
c             v(f0)) is used in testing for relative function
c             convergence.
c v(reldx).... v(17) is the scaled relative change in x caused by the
c             current step, computed as
c                  max(abs(d(i)*(x(i)-x0(i)), 1 .le. i .le. p) /
c                     max(d(i)*(abs(x(i))+abs(x0(i))), 1 .le. i .le. p),
c             where x = x0 + step.
c
c-------------------------------  notes  -------------------------------
c
c  ***  algorithm notes  ***
c
c        see ref. 1 for a description of the algorithm used.
c        on problems which are naturally well scaled, better perform-
c     ance may be obtained by setting v(d0init) = 1.0 and iv(dtype) = 0,
c     which will cause the scale vector d to be set to all ones.
c
c  ***  usage notes  ***
c
c        after a return with iv(1) .le. 11, it is possible to restart,
c     i.e., to change some of the iv and v input values described above
c     and continue the algorithm from the point where it was interrupt-
c     ed.  iv(1) should not be changed, nor should any entries of iv
c     and v other than the input values (those supplied by dfault).
c        those who do not wish to write a calcj which computes the ja-
c     cobian matrix analytically should call nl2sno rather than nl2sol.
c     nl2sno uses finite differences to compute an approximate jacobian.
c        those who would prefer to provide r and j (the residual and
c     jacobian) by reverse communication rather than by writing subrou-
c     tines calcr and calcj may call on nl2itr directly.  see the com-
c     ments at the beginning of nl2itr.
c        those who use nl2sol interactively may wish to supply their
c     own stopx function, which should return .true. if the break key
c     has been pressed since stopx was last invoked.  this makes it pos-
c     sible to externally interrupt nl2sol (which will return with
c     iv(1) = 11 if stopx returns .true.).
c        storage for j is allocated at the end of v.  thus the caller
c     may make v longer than specified above and may allow calcj to use
c     elements of j beyond the first n*p as scratch storage.
c
c  ***  portability notes  ***
c
c        the nl2sol distribution tape contains both single- and double-
c     precision versions of the nl2sol source code, so it should be un-
c     necessary to change precisions.
c        only the functions imdcon and rmdcon contain machine-dependent
c     constants.  to change from one machine to another, it should
c     suffice to change the (few) relevant lines in these functions.
c        intrinsic functions are explicitly declared.  on certain com-
c     puters (e.g. univac), it may be necessary to comment out these
c     declarations.  so that this may be done automatically by a simple
c     program, such declarations are preceded by a comment having c/+
c     in columns 1-3 and blanks in columns 4-72 and are followed by
c     a comment having c/ in columns 1 and 2 and blanks in columns 3-72.
c        the nl2sol source code is expressed in 1966 ansi standard
c     fortran.  it may be converted to fortran 77 by
c     commenting out all lines that fall between a line having c/6 in
c     columns 1-3 and a line having c/7 in columns 1-3 and by removing
c     (i.e., replacing by a blank) the c in column 1 of the lines that
c     follow the c/7 line and preceed a line having c/ in columns 1-2
c     and blanks in columns 3-72.  these changes convert some data
c     statements into parameter statements, convert some variables from
c     real to character*4, and make the data statements that initialize
c     these variables use character strings delimited by primes instead
c     of hollerith constants.  (such variables and data statements
c     appear only in modules itsmry and parchk.  parameter statements
c     appear nearly everywhere.)
c
c  ***  references  ***
c
c 1.  dennis, j.e., gay, d.m., and welsch, r.e. (1981), an adaptive
c             nonlinear least-squares algorithm, acm trans. math.
c             software, vol. 7, no. 3.
c
c
c  ***  general  ***
c
c     coded by david m. gay (winter 1979 - winter 1980).
c     this subroutine was written in connection with research
c     supported by the national science foundation under grants
c     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and
c     mcs-7906671.
c
c----------------------------  declarations  ---------------------------
c
      external itsmry, nl2itr
c itsmry... prints iteration summary and info about initial and final x.
c nl2itr... reverse-communication routine that carries out nl2sol algo-
c             rithm.
c
      logical strted
      integer d1, j1, nf, r1
c
c  ***  subscripts for iv and v  ***
c
      integer d, j, nfcall, nfgcal, r, toobig
c
c  ***  iv subscript values  ***
c
c/6
      data nfcall/6/, nfgcal/7/, toobig/2/
c/7
c     parameter (nfcall=6, nfgcal=7, toobig=2)
c/
c
c  ***  v subscript values  ***
c
c/6
      data d/27/, j/33/, r/50/
c/7
c     parameter (d=27, j=33, r=50)
c/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
      d1 = 94 + 2*n + p*(3*p + 31)/2
      iv(d) = d1
      r1 = d1 + p
      iv(r) = r1
      j1 = r1 + n
      iv(j) = j1
      strted = .true.
      if (iv(1) .ne. 0 .and. iv(1) .ne. 12) go to 40
         strted = .false.
         iv(nfcall) = 1
         iv(nfgcal) = 1
c
 10   nf = iv(nfcall)
      call calcr(n, p, x, nf, v(r1), uiparm, urparm, ufparm)
      if (strted) go to 20
         if (nf .gt. 0) go to 30
              iv(1) = 13
              go to 60
c
 20   if (nf .le. 0) iv(toobig) = 1
      go to 40
c
 30   call calcj(n, p, x, iv(nfgcal), v(j1), uiparm, urparm, ufparm)
      if (iv(nfgcal) .eq. 0) go to 50
      strted = .true.
c
 40   call nl2itr(v(d1), iv, v(j1), n, n, p, v(r1), v, x)
      if (iv(1) - 2) 10, 30, 999
c
 50   iv(1) = 15
 60   call itsmry(v(d1), iv, p, v, x)
c
 999  return
c  ***  last card of nl2sol follows  ***
      end
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
      data hfac/1.d+3/, negpt5/-0.5d+0/, one/1.d+0/, zero/0.d+0/
c/7
c     parameter (hfac=1.d+3, negpt5=-0.5d+0, one=1.d+0, zero=0.d+0)
c/
c
c  ***  iv subscript values  ***
c
c/6
      data covprt/14/, covreq/15/, d/27/, dtype/16/, j/33/,
     1     nfcall/6/, nfgcal/7/, r/50/, toobig/2/
c/7
c     parameter (covprt=14, covreq=15, d=27, dtype=16, j=33,
c    1     nfcall=6, nfgcal=7, r=50, toobig=2)
c/
c
c  ***  v subscript values  ***
c
c/6
      data dltfdj/36/, dinit/38/
c/7
c     parameter (dltfdj=36)
c     save hlim
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
       if (v(dinit).gt.zero) call vscopy(p, v(d1), v(dinit))
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
      data cnvcod/34/, covmat/26/, covprt/14/,
     1     covreq/15/, dig/43/, dtype/16/, g/28/, h/44/,
     2     ierr/32/, inits/25/, ipivot/61/, ipiv0/60/,
     3     irc/3/, kagqt/35/, kalm/36/, lky/37/, lmat/58/,
     4     mode/38/, model/5/, mxfcal/17/, mxiter/18/,
     5     nfcall/6/, nfgcal/7/, nfcov/40/, ngcov/41/,
     6     ngcall/30/, niter/31/, qtr/49/,
     7     radinc/8/, rd/51/, restor/9/, rsave/52/, s/53/,
     8     step/55/, stglim/11/, stlstg/56/, sused/57/,
     9     switch/12/, toobig/2/, w/59/, xirc/13/, x0/60/
c/7
c     parameter (cnvcod=34, covmat=26, covprt=14,
c    1     covreq=15, dig=43, dtype=16, g=28, h=44,
c    2     ierr=32, inits=25, ipivot=61, ipiv0=60,
c    3     irc=3, kagqt=35, kalm=36, lky=37, lmat=58,
c    4     mode=38, model=5, mxfcal=17, mxiter=18,
c    5     nfcall=6, nfgcal=7, nfcov=40, ngcov=41,
c    6     ngcall=30, niter=31, qtr=49,
c    7     radinc=8, rd=51, restor=9, rsave=52, s=53,
c    8     step=55, stglim=11, stlstg=56, sused=57,
c    9     switch=12, toobig=2, w=59, xirc=13, x0=60)
c/
c
c  ***  v subscript values  ***
c
c/6
      data cosmin/43/, dgnorm/1/, dinit/38/, dstnrm/2/,
     1     d0init/37/, f/10/, fdif/11/, fuzz/45/,
     2     f0/13/, gtstep/4/, incfac/23/,
     3     jtinit/39/, jtol1/87/, lmax0/35/,
     4     nvsave/9/, phmxfc/21/, preduc/7/,
     5     radfac/16/, radius/8/, rad0/9/, rlimit/42/,
     6     size/47/, stppar/5/, tuner4/29/, tuner5/30/,
     7     vsave1/78/, wscale/48/
c/7
c     parameter (cosmin=43, dgnorm=1, dinit=38, dstnrm=2,
c    1     d0init=37, f=10, fdif=11, fuzz=45,
c    2     f0=13, gtstep=4, incfac=23,
c    3     jtinit=39, jtol1=87, lmax0=35,
c    4     nvsave=9, phmxfc=21, preduc=7,
c    5     radfac=16, radius=8, rad0=9, rlimit=42,
c    6     size=47, stppar=5, tuner4=29, tuner5=30,
c    7     vsave1=78, wscale=48)
c/
c
c
c/6
      data half/0.5d+0/, negone/-1.d+0/, one/1.d+0/, zero/0.d+0/
c/7
c     parameter (half=0.5d+0, negone=-1.d+0, one=1.d+0, zero=0.d+0)
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
      subroutine assess (d, iv, p, step, stlstg, v, x, x0)
c
c  ***  assess candidate step (nl2sol version 2.2)  ***
c
      integer p, iv(13)
      double precision d(p), step(p), stlstg(p), v(35), x(p), x0(p)
c
c  ***  purpose  ***
c
c        this subroutine is called by an unconstrained minimization
c     routine to assess the next candidate step.  it may recommend one
c     of several courses of action, such as accepting the step, recom-
c     puting it using the same or a new quadratic model, or halting due
c     to convergence or false convergence.  see the return code listing
c     below.
c
c--------------------------  parameter usage  --------------------------
c
c     iv (i/o) integer parameter and scratch vector -- see description
c             below of iv values referenced.
c      d (in)  scale vector used in computing v(reldx) -- see below.
c      p (in)  number of parameters being optimized.
c   step (i/o) on input, step is the step to be assessed.  it is un-
c             changed on output unless a previous step achieved a
c             better objective function reduction, in which case stlstg
c             will have been copied to step.
c stlstg (i/o) when assess recommends recomputing step even though the
c             current (or a previous) step yields an objective func-
c             tion decrease, it saves in stlstg the step that gave the
c             best function reduction seen so far (in the current itera-
c             tion).  if the recomputed step yields a larger function
c             value, then step is restored from stlstg and
c             x = x0 + step is recomputed.
c      v (i/o) real parameter and scratch vector -- see description
c             below of v values referenced.
c      x (i/o) on input, x = x0 + step is the point at which the objec-
c             tive function has just been evaluated.  if an earlier
c             step yielded a bigger function decrease, then x is
c             restored to the corresponding earlier value.  otherwise,
c             if the current step does not give any function decrease,
c             then x is restored to x0.
c     x0 (in)  initial objective function parameter vector (at the
c             start of the current iteration).
c
c  ***  iv values referenced  ***
c
c    iv(irc) (i/o) on input for the first step tried in a new iteration,
c             iv(irc) should be set to 3 or 4 (the value to which it is
c             set when step is definitely to be accepted).  on input
c             after step has been recomputed, iv(irc) should be
c             unchanged since the previous return of assess.
c                on output, iv(irc) is a return code having one of the
c             following values...
c                  1 = switch models or try smaller step.
c                  2 = switch models or accept step.
c                  3 = accept step and determine v(radfac) by gradient
c                       tests.
c                  4 = accept step, v(radfac) has been determined.
c                  5 = recompute step (using the same model).
c                  6 = recompute step with radius = v(lmax0) but do not
c                       evaulate the objective function.
c                  7 = x-convergence (see v(xctol)).
c                  8 = relative function convergence (see v(rfctol)).
c                  9 = both x- and relative function convergence.
c                 10 = absolute function convergence (see v(afctol)).
c                 11 = singular convergence (see v(lmax0)).
c                 12 = false convergence (see v(xftol)).
c                 13 = iv(irc) was out of range on input.
c             return code i has precdence over i+1 for i = 9, 10, 11.
c iv(mlstgd) (i/o) saved value of iv(model).
c  iv(model) (i/o) on input, iv(model) should be an integer identifying
c             the current quadratic model of the objective function.
c             if a previous step yielded a better function reduction,
c             then iv(model) will be set to iv(mlstgd) on output.
c iv(nfcall) (in)  invocation count for the objective function.
c iv(nfgcal) (i/o) value of iv(nfcall) at step that gave the biggest
c             function reduction this iteration.  iv(nfgcal) remains
c             unchanged until a function reduction is obtained.
c iv(radinc) (i/o) the number of radius increases (or minus the number
c             of decreases) so far this iteration.
c iv(restor) (out) set to 0 unless x and v(f) have been restored, in
c             which case assess sets iv(restor) = 1.
c  iv(stage) (i/o) count of the number of models tried so far in the
c             current iteration.
c iv(stglim) (in)  maximum number of models to consider.
c iv(switch) (out) set to 0 unless a new model is being tried and it
c             gives a smaller function value than the previous model,
c             in which case assess sets iv(switch) = 1.
c iv(toobig) (in)  is nonzero if step was too big (e.g. if it caused
c             overflow).
c   iv(xirc) (i/o) value that iv(irc) would have in the absence of
c             convergence, false convergence, and oversized steps.
c
c  ***  v values referenced  ***
c
c v(afctol) (in)  absolute function convergence tolerance.  if the
c             absolute value of the current function value v(f) is less
c             than v(afctol), then assess returns with iv(irc) = 10.
c v(decfac) (in)  factor by which to decrease radius when iv(toobig) is
c             nonzero.
c v(dstnrm) (in)  the 2-norm of d*step.
c v(dstsav) (i/o) value of v(dstnrm) on saved step.
c   v(dst0) (in)  the 2-norm of d times the newton step (when defined,
c             i.e., for v(nreduc) .ge. 0).
c      v(f) (i/o) on both input and output, v(f) is the objective func-
c             tion value at x.  if x is restored to a previous value,
c             then v(f) is restored to the corresponding value.
c   v(fdif) (out) the function reduction v(f0) - v(f) (for the output
c             value of v(f) if an earlier step gave a bigger function
c             decrease, and for the input value of v(f) otherwise).
c v(flstgd) (i/o) saved value of v(f).
c     v(f0) (in)  objective function value at start of iteration.
c v(gtslst) (i/o) value of v(gtstep) on saved step.
c v(gtstep) (in)  inner product between step and gradient.
c v(incfac) (in)  minimum factor by which to increase radius.
c  v(lmax0) (in)  maximum reasonable step size (and initial step bound).
c             if the actual function decrease is no more than twice
c             what was predicted, if a return with iv(irc) = 7, 8, 9,
c             or 10 does not occur, if v(dstnrm) .gt. v(lmax0), and if
c             v(preduc) .le. v(rfctol) * abs(v(f0)), then assess re-
c             turns with iv(irc) = 11.  if so doing appears worthwhile,
c             then assess repeats this test with v(preduc) computed for
c             a step of length v(lmax0) (by a return with iv(irc) = 6).
c v(nreduc) (i/o)  function reduction predicted by quadratic model for
c             newton step.  if assess is called with iv(irc) = 6, i.e.,
c             if v(preduc) has been computed with radius = v(lmax0) for
c             use in the singular convervence test, then v(nreduc) is
c             set to -v(preduc) before the latter is restored.
c v(plstgd) (i/o) value of v(preduc) on saved step.
c v(preduc) (i/o) function reduction predicted by quadratic model for
c             current step.
c v(radfac) (out) factor to be used in determining the new radius,
c             which should be v(radfac)*dst, where  dst  is either the
c             output value of v(dstnrm) or the 2-norm of
c             diag(newd)*step  for the output value of step and the
c             updated version, newd, of the scale vector d.  for
c             iv(irc) = 3, v(radfac) = 1.0 is returned.
c v(rdfcmn) (in)  minimum value for v(radfac) in terms of the input
c             value of v(dstnrm) -- suggested value = 0.1.
c v(rdfcmx) (in)  maximum value for v(radfac) -- suggested value = 4.0.
c  v(reldx) (out) scaled relative change in x caused by step, computed
c             by function  reldst  as
c                 max (d(i)*abs(x(i)-x0(i)), 1 .le. i .le. p) /
c                    max (d(i)*(abs(x(i))+abs(x0(i))), 1 .le. i .le. p).
c             if an acceptable step is returned, then v(reldx) is com-
c             puted using the output (possibly restored) values of x
c             and step.  otherwise it is computed using the input
c             values.
c v(rfctol) (in)  relative function convergence tolerance.  if the
c             actual function reduction is at most twice what was pre-
c             dicted and  v(nreduc) .le. v(rfctol)*abs(v(f0)),  then
c             assess returns with iv(irc) = 8 or 9.  see also v(lmax0).
c v(stppar) (in)  marquardt parameter -- 0 means full newton step.
c v(tuner1) (in)  tuning constant used to decide if the function
c             reduction was much less than expected.  suggested
c             value = 0.1.
c v(tuner2) (in)  tuning constant used to decide if the function
c             reduction was large enough to accept step.  suggested
c             value = 10**-4.
c v(tuner3) (in)  tuning constant used to decide if the radius
c             should be increased.  suggested value = 0.75.
c  v(xctol) (in)  x-convergence criterion.  if step is a newton step
c             (v(stppar) = 0) having v(reldx) .le. v(xctol) and giving
c             at most twice the predicted function decrease, then
c             assess returns iv(irc) = 7 or 9.
c  v(xftol) (in)  false convergence tolerance.  if step gave no or only
c             a small function decrease and v(reldx) .le. v(xftol),
c             then assess returns with iv(irc) = 12.
c
c-------------------------------  notes  -------------------------------
c
c  ***  application and usage restrictions  ***
c
c        this routine is called as part of the nl2sol (nonlinear
c     least-squares) package.  it may be used in any unconstrained
c     minimization solver that uses dogleg, goldfeld-quandt-trotter,
c     or levenberg-marquardt steps.
c
c  ***  algorithm notes  ***
c
c        see (1) for further discussion of the assessing and model
c     switching strategies.  while nl2sol considers only two models,
c     assess is designed to handle any number of models.
c
c  ***  usage notes  ***
c
c        on the first call of an iteration, only the i/o variables
c     step, x, iv(irc), iv(model), v(f), v(dstnrm), v(gtstep), and
c     v(preduc) need have been initialized.  between calls, no i/o
c     values execpt step, x, iv(model), v(f) and the stopping toler-
c     ances should be changed.
c        after a return for convergence or false convergence, one can
c     change the stopping tolerances and call assess again, in which
c     case the stopping tests will be repeated.
c
c  ***  references  ***
c
c     (1) dennis, j.e., jr., gay, d.m., and welsch, r.e. (1981),
c        an adaptive nonlinear least-squares algorithm,
c        acm trans. math. software, vol. 7, no. 3.
c
c     (2) powell, m.j.d. (1970)  a fortran subroutine for solving
c        systems of nonlinear algebraic equations, in numerical
c        methods for nonlinear algebraic equations, edited by
c        p. rabinowitz, gordon and breach, london.
c
c  ***  history  ***
c
c        john dennis designed much of this routine, starting with
c     ideas in (2). roy welsch suggested the model switching strategy.
c        david gay and stephen peters cast this subroutine into a more
c     portable form (winter 1977), and david gay cast it into its
c     present form (fall 1978).
c
c  ***  general  ***
c
c     this subroutine was written in connection with research
c     supported by the national science foundation under grants
c     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and
c     mcs-7906671.
c
c------------------------  external quantities  ------------------------
c
c  ***  external functions and subroutines  ***
c
      external reldst, vcopy
      double precision reldst
c
c vcopy.... copies one vector to another.
c
c  ***  intrinsic functions  ***
c/+
      integer iabs
      double precision dabs, dmax1
c/
c  ***  no common blocks  ***
c
c--------------------------  local variables  --------------------------
c
      logical goodx
      integer i, nfc
      double precision emax, gts, half, one, reldx1, rfac1, two, xmax,
     1                 zero
c
c  ***  subscripts for iv and v  ***
c
      integer afctol, decfac, dstnrm, dstsav, dst0, f, fdif, flstgd, f0,
     1        gtslst, gtstep, incfac, irc, lmax0, mlstgd, model, nfcall,
     2        nfgcal, nreduc, plstgd, preduc, radfac, radinc, rdfcmn,
     3        rdfcmx, reldx, restor, rfctol, stage, stglim, stppar,
     4        switch, toobig, tuner1, tuner2, tuner3, xctol, xftol,
     5        xirc
c
c  ***  data initializations  ***
c
c/6
      data half/0.5d+0/, one/1.d+0/, two/2.d+0/, zero/0.d+0/
c/7
c     parameter (half=0.5d+0, one=1.d+0, two=2.d+0, zero=0.d+0)
c/
c
c/6
      data irc/3/, mlstgd/4/, model/5/, nfcall/6/,
     1     nfgcal/7/, radinc/8/, restor/9/, stage/10/,
     2     stglim/11/, switch/12/, toobig/2/, xirc/13/
c/7
c     parameter (irc=3, mlstgd=4, model=5, nfcall=6,
c    1     nfgcal=7, radinc=8, restor=9, stage=10,
c    2     stglim=11, switch=12, toobig=2, xirc=13)
c/
c/6
      data afctol/31/, decfac/22/, dstnrm/2/, dst0/3/,
     1     dstsav/18/, f/10/, fdif/11/, flstgd/12/, f0/13/,
     2     gtslst/14/, gtstep/4/, incfac/23/,
     3     lmax0/35/, nreduc/6/, plstgd/15/, preduc/7/,
     4     radfac/16/, rdfcmn/24/, rdfcmx/25/,
     5     reldx/17/, rfctol/32/, stppar/5/, tuner1/26/,
     6     tuner2/27/, tuner3/28/, xctol/33/, xftol/34/
c/7
c     parameter (afctol=31, decfac=22, dstnrm=2, dst0=3,
c    1     dstsav=18, f=10, fdif=11, flstgd=12, f0=13,
c    2     gtslst=14, gtstep=4, incfac=23,
c    3     lmax0=35, nreduc=6, plstgd=15, preduc=7,
c    4     radfac=16, rdfcmn=24, rdfcmx=25,
c    5     reldx=17, rfctol=32, stppar=5, tuner1=26,
c    6     tuner2=27, tuner3=28, xctol=33, xftol=34)
c/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
      nfc = iv(nfcall)
      iv(switch) = 0
      iv(restor) = 0
      rfac1 = one
      goodx = .true.
      i = iv(irc)
      if (i .ge. 1 .and. i .le. 12)
     1             go to (20,30,10,10,40,360,290,290,290,290,290,140), i
         iv(irc) = 13
         go to 999
c
c  ***  initialize for new iteration  ***
c
 10   iv(stage) = 1
      iv(radinc) = 0
      v(flstgd) = v(f0)
      if (iv(toobig) .eq. 0) go to 90
         iv(stage) = -1
         iv(xirc) = i
         go to 60
c
c  ***  step was recomputed with new model or smaller radius  ***
c  ***  first decide which  ***
c
 20   if (iv(model) .ne. iv(mlstgd)) go to 30
c        ***  old model retained, smaller radius tried  ***
c        ***  do not consider any more new models this iteration  ***
         iv(stage) = iv(stglim)
         iv(radinc) = -1
         go to 90
c
c  ***  a new model is being tried.  decide whether to keep it.  ***
c
 30   iv(stage) = iv(stage) + 1
c
c     ***  now we add the possibiltiy that step was recomputed with  ***
c     ***  the same model, perhaps because of an oversized step.     ***
c
 40   if (iv(stage) .gt. 0) go to 50
c
c        ***  step was recomputed because it was too big.  ***
c
         if (iv(toobig) .ne. 0) go to 60
c
c        ***  restore iv(stage) and pick up where we left off.  ***
c
         iv(stage) = -iv(stage)
         i = iv(xirc)
         go to (20, 30, 90, 90, 70), i
c
 50   if (iv(toobig) .eq. 0) go to 70
c
c  ***  handle oversize step  ***
c
      if (iv(radinc) .gt. 0) go to 80
         iv(stage) = -iv(stage)
         iv(xirc) = iv(irc)
c
 60      v(radfac) = v(decfac)
         iv(radinc) = iv(radinc) - 1
         iv(irc) = 5
         go to 999
c
 70   if (v(f) .lt. v(flstgd)) go to 90
c
c     *** the new step is a loser.  restore old model.  ***
c
      if (iv(model) .eq. iv(mlstgd)) go to 80
         iv(model) = iv(mlstgd)
         iv(switch) = 1
c
c     ***  restore step, etc. only if a previous step decreased v(f).
c
 80   if (v(flstgd) .ge. v(f0)) go to 90
         iv(restor) = 1
         v(f) = v(flstgd)
         v(preduc) = v(plstgd)
         v(gtstep) = v(gtslst)
         if (iv(switch) .eq. 0) rfac1 = v(dstnrm) / v(dstsav)
         v(dstnrm) = v(dstsav)
         nfc = iv(nfgcal)
         goodx = .false.
c
c
c  ***  compute relative change in x by current step  ***
c
 90   reldx1 = reldst(p, d, x, x0)
c
c  ***  restore x and step if necessary  ***
c
      if (goodx) go to 105
      do 100 i = 1, p
         step(i) = stlstg(i)
         x(i) = x0(i) + stlstg(i)
 100     continue
c
 105  v(fdif) = v(f0) - v(f)
      if (v(fdif) .gt. v(tuner2) * v(preduc)) go to 120
c
c        ***  no (or only a trivial) function decrease
c        ***  -- so try new model or smaller radius
c
         v(reldx) = reldx1
         if (v(f) .lt. v(f0)) go to 110
              iv(mlstgd) = iv(model)
              v(flstgd) = v(f)
              v(f) = v(f0)
              call vcopy(p, x, x0)
              iv(restor) = 1
              go to 115
 110     iv(nfgcal) = nfc
 115     iv(irc) = 1
         if (iv(stage) .lt. iv(stglim)) go to 130
              iv(irc) = 5
              iv(radinc) = iv(radinc) - 1
              go to 130
c
c  ***  nontrivial function decrease achieved  ***
c
 120  iv(nfgcal) = nfc
      rfac1 = one
      if (goodx) v(reldx) = reldx1
      v(dstsav) = v(dstnrm)
      if (v(fdif) .gt. v(preduc)*v(tuner1)) go to 200
c
c  ***  decrease was much less than predicted -- either change models
c  ***  or accept step with decreased radius.
c
      if (iv(stage) .ge. iv(stglim)) go to 125
c        ***  consider switching models  ***
         iv(irc) = 2
         go to 130
c
c     ***  accept step with decreased radius  ***
c
 125  iv(irc) = 4
c
c  ***  set v(radfac) to fletcher*s decrease factor  ***
c
 130  iv(xirc) = iv(irc)
      emax = v(gtstep) + v(fdif)
      v(radfac) = half * rfac1
      if (emax .lt. v(gtstep)) v(radfac) = rfac1 * dmax1(v(rdfcmn),
     1                                           half * v(gtstep)/emax)
c
c  ***  do false convergence test  ***
c
 140  if (v(reldx) .le. v(xftol)) go to 160
         iv(irc) = iv(xirc)
         if (v(f) .lt. v(f0)) go to 230
              go to 300
c
 160  iv(irc) = 12
      go to 310
c
c  ***  handle good function decrease  ***
c
 200  if (v(fdif) .lt. (-v(tuner3) * v(gtstep))) go to 260
c
c     ***  increasing radius looks worthwhile.  see if we just
c     ***  recomputed step with a decreased radius or restored step
c     ***  after recomputing it with a larger radius.
c
      if (iv(radinc) .lt. 0) go to 260
      if (iv(restor) .eq. 1) go to 260
c
c        ***  we did not.  try a longer step unless this was a newton
c        ***  step.
c
         v(radfac) = v(rdfcmx)
         gts = v(gtstep)
         if (v(fdif) .lt. (half/v(radfac) - one) * gts)
     1            v(radfac) = dmax1(v(incfac), half*gts/(gts + v(fdif)))
         iv(irc) = 4
         if (v(stppar) .eq. zero) go to 300
c             ***  step was not a newton step.  recompute it with
c             ***  a larger radius.
              iv(irc) = 5
              iv(radinc) = iv(radinc) + 1
c
c  ***  save values corresponding to good step  ***
c
 230  v(flstgd) = v(f)
      iv(mlstgd) = iv(model)
      call vcopy(p, stlstg, step)
      v(dstsav) = v(dstnrm)
      iv(nfgcal) = nfc
      v(plstgd) = v(preduc)
      v(gtslst) = v(gtstep)
      go to 300
c
c  ***  accept step with radius unchanged  ***
c
 260  v(radfac) = one
      iv(irc) = 3
      go to 300
c
c  ***  come here for a restart after convergence  ***
c
 290  iv(irc) = iv(xirc)
      if (v(dstsav) .ge. zero) go to 310
         iv(irc) = 12
         go to 310
c
c  ***  perform convergence tests  ***
c
 300  iv(xirc) = iv(irc)
 310  if (dabs(v(f)) .lt. v(afctol)) iv(irc) = 10
      if (half * v(fdif) .gt. v(preduc)) go to 999
      emax = v(rfctol) * dabs(v(f0))
      if (v(dstnrm) .gt. v(lmax0) .and. v(preduc) .le. emax)
     1                       iv(irc) = 11
      if (v(dst0) .lt. zero) go to 320
      i = 0
      if ((v(nreduc) .gt. zero .and. v(nreduc) .le. emax) .or.
     1    (v(nreduc) .eq. zero. and. v(preduc) .eq. zero))  i = 2
      if (v(stppar) .eq. zero .and. v(reldx) .le. v(xctol)
     1                        .and. goodx)                  i = i + 1
      if (i .gt. 0) iv(irc) = i + 6
c
c  ***  consider recomputing step of length v(lmax0) for singular
c  ***  convergence test.
c
 320  if (iabs(iv(irc)-3) .gt. 2 .and. iv(irc) .ne. 12) go to 999
      if (v(dstnrm) .gt. v(lmax0)) go to 330
         if (v(preduc) .ge. emax) go to 999
              if (v(dst0) .le. zero) go to 340
                   if (half * v(dst0) .le. v(lmax0)) go to 999
                        go to 340
 330  if (half * v(dstnrm) .le. v(lmax0)) go to 999
      xmax = v(lmax0) / v(dstnrm)
      if (xmax * (two - xmax) * v(preduc) .ge. emax) go to 999
 340  if (v(nreduc) .lt. zero) go to 370
c
c  ***  recompute v(preduc) for use in singular convergence test  ***
c
      v(gtslst) = v(gtstep)
      v(dstsav) = v(dstnrm)
      if (iv(irc) .eq. 12) v(dstsav) = -v(dstsav)
      v(plstgd) = v(preduc)
      iv(irc) = 6
      call vcopy(p, stlstg, step)
      go to 999
c
c  ***  perform singular convergence test with recomputed v(preduc)  ***
c
 360  v(gtstep) = v(gtslst)
      v(dstnrm) = dabs(v(dstsav))
      call vcopy(p, step, stlstg)
      iv(irc) = iv(xirc)
      if (v(dstsav) .le. zero) iv(irc) = 12
      v(nreduc) = -v(preduc)
      v(preduc) = v(plstgd)
 370  if (-v(nreduc) .le. v(rfctol) * dabs(v(f0))) iv(irc) = 11
c
 999  return
c
c  ***  last card of assess follows  ***
      end
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
      data half/0.5d+0/, negpt5/-0.5d+0/, one/1.d+0/, two/2.d+0/,
     1     zero/0.d+0/
c/7
c     parameter (half=0.5d+0, negpt5=-0.5d+0, one=1.d+0, two=2.d+0,
c    1     zero=0.d+0)
c/
c
c/6
      data covmat/26/, covreq/15/, delta/50/, delta0/44/,
     1     dltfdc/40/, f/10/, fx/46/, g/28/, h/44/, ierr/32/,
     2     ipivot/61/, ipiv0/60/, kagqt/35/, kalm/36/,
     3     lmat/58/, mode/38/, nfgcal/7/, qtr/49/,
     4     rd/51/, rsave/52/, savei/54/, switch/12/,
     5     toobig/2/, w/59/, xmsave/49/
c/7
c     parameter (covmat=26, covreq=15, delta=50, delta0=44,
c    1     dltfdc=40, f=10, fx=46, g=28, h=44, ierr=32,
c    2     ipivot=61, ipiv0=60, kagqt=35, kalm=36,
c    3     lmat=58, mode=38, nfgcal=7, qtr=49,
c    4     rd=51, rsave=52, savei=54, switch=12,
c    5     toobig=2, w=59, xmsave=49)
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
      subroutine dfault(iv, v)
c
c  ***  supply nl2sol (version 2.2) default values to iv and v  ***
c
      integer iv(25)
      double precision v(45)
c/+
      double precision dmax1
c/
      external imdcon, rmdcon
      integer imdcon
      double precision rmdcon
c
      double precision machep, mepcrt, one, sqteps, three
c
c  ***  subscripts for iv and v  ***
c
      integer afctol, cosmin, covprt, covreq, decfac, delta0, dfac,
     1        dinit, dltfdc, dltfdj, dtype, d0init, epslon, fuzz,
     2        incfac, inits, jtinit, lmax0, mxfcal, mxiter, outlev,
     3        parprt, phmnfc, phmxfc, prunit, rdfcmn, rdfcmx,
     4        rfctol, rlimit, solprt, statpr, tuner1, tuner2, tuner3,
     5        tuner4, tuner5, xctol, xftol, x0prt
c
c/6
      data one/1.d+0/, three/3.d+0/
c/7
c     parameter (one=1.d+0, three=3.d+0)
c/
c
c  ***  iv subscript values  ***
c
c/6
      data covprt/14/, covreq/15/, dtype/16/, inits/25/,
     1     mxfcal/17/, mxiter/18/, outlev/19/,
     2     parprt/20/, prunit/21/, solprt/22/,
     3     statpr/23/, x0prt/24/
c/7
c     parameter (covprt=14, covreq=15, dtype=16, inits=25,
c    1     mxfcal=17, mxiter=18, outlev=19,
c    2     parprt=20, prunit=21, solprt=22,
c    3     statpr=23, x0prt=24)
c/
c
c  ***  v subscript values  ***
c
c/6
      data afctol/31/, cosmin/43/, decfac/22/,
     1     delta0/44/, dfac/41/, dinit/38/, dltfdc/40/,
     2     dltfdj/36/, d0init/37/, epslon/19/, fuzz/45/,
     3     incfac/23/, jtinit/39/, lmax0/35/, phmnfc/20/,
     4     phmxfc/21/, rdfcmn/24/, rdfcmx/25/,
     5     rfctol/32/, rlimit/42/, tuner1/26/,
     6     tuner2/27/, tuner3/28/, tuner4/29/,
     7     tuner5/30/, xctol/33/, xftol/34/
c/7
c     parameter (afctol=31, cosmin=43, decfac=22,
c    1     delta0=44, dfac=41, dinit=38, dltfdc=40,
c    2     dltfdj=36, d0init=37, epslon=19, fuzz=45,
c    3     incfac=23, jtinit=39, lmax0=35, phmnfc=20,
c    4     phmxfc=21, rdfcmn=24, rdfcmx=25,
c    5     rfctol=32, rlimit=42, tuner1=26,
c    6     tuner2=27, tuner3=28, tuner4=29,
c    7     tuner5=30, xctol=33, xftol=34)
c/
c
c-----------------------------------------------------------------------
c
      iv(1) = 12
      iv(covprt) = 1
      iv(covreq) = 1
      iv(dtype) = 1
      iv(inits) = 0
      iv(mxfcal) = 200
      iv(mxiter) = 150
      iv(outlev) = 1
      iv(parprt) = 1
      iv(prunit) = imdcon(1)
      iv(solprt) = 1
      iv(statpr) = 1
      iv(x0prt) = 1
c
      machep = rmdcon(3)
      v(afctol) = 1.d-20
      if (machep .gt. 1.d-10) v(afctol) = machep**2
      v(cosmin) = dmax1(1.d-6, 1.d+2 * machep)
      v(decfac) = 0.5d+0
      sqteps = rmdcon(4)
      v(delta0) = sqteps
      v(dfac) = 0.6d+0
      v(dinit) = 0.d+0
      mepcrt = machep ** (one/three)
      v(dltfdc) = mepcrt
      v(dltfdj) = sqteps
      v(d0init) = 1.d+0
      v(epslon) = 0.1d+0
      v(fuzz) = 1.5d+0
      v(incfac) = 2.d+0
      v(jtinit) = 1.d-6
      v(lmax0) = 100.d+0
      v(phmnfc) = -0.1d+0
      v(phmxfc) = 0.1d+0
      v(rdfcmn) = 0.1d+0
      v(rdfcmx) = 4.d+0
      v(rfctol) = dmax1(1.d-10, mepcrt**2)
      v(rlimit) = rmdcon(5)
      v(tuner1) = 0.1d+0
      v(tuner2) = 1.d-4
      v(tuner3) = 0.75d+0
      v(tuner4) = 0.5d+0
      v(tuner5) = 0.75d+0
      v(xctol) = sqteps
      v(xftol) = 1.d+2 * machep
c
 999  return
c  ***  last card of dfault follows  ***
      end
      double precision function dotprd(p, x, y)
c
c  ***  return the inner product of the p-vectors x and y.  ***
c
      integer p
      double precision x(p), y(p)
c
      integer i
      double precision one, sqteta, t, zero
c/+
      double precision dmax1, dabs
c/
      external rmdcon
      double precision rmdcon
c
c  ***  rmdcon(2) returns a machine-dependent constant, sqteta, which
c  ***  is slightly larger than the smallest positive number that
c  ***  can be squared without underflowing.
c
c/6
      data one/1.d+0/, sqteta/0.d+0/, zero/0.d+0/
c/7
c     parameter (one=1.d+0, zero=0.d+0)
c     data sqteta/0.d+0/
c/
c
      dotprd = zero
      if (p .le. 0) go to 999
      if (sqteta .eq. zero) sqteta = rmdcon(2)
      do 20 i = 1, p
         t = dmax1(dabs(x(i)), dabs(y(i)))
         if (t .gt. one) go to 10
         if (t .lt. sqteta) go to 20
         t = (x(i)/sqteta)*y(i)
         if (dabs(t) .lt. sqteta) go to 20
 10      dotprd = dotprd + x(i)*y(i)
 20   continue
c
 999  return
c  ***  last card of dotprd follows  ***
      end
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
      data dfac/41/, dtype/16/, jtol0/86/, niter/31/, s/53/
c/7
c     parameter (dfac=41, dtype=16, jtol0=86, niter=31, s=53)
c/
c
c/6
      data zero/0.d+0/
c/7
c     parameter (zero=0.d+0)
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
      subroutine gqtstp(d, dig, dihdi, ka, l, p, step, v, w)
c
c  *** compute goldfeld-quandt-trotter step by more-hebden technique ***
c  ***  (nl2sol version 2.2)  ***
c
c  ***  parameter declarations  ***
c
      integer ka, p
      double precision d(p), dig(p), dihdi(1), l(1), v(21), step(p),
     1                 w(1)
c     dimension dihdi(p*(p+1)/2), l(p*(p+1)/2), w(4*p+7)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  purpose  ***
c
c        given the (compactly stored) lower triangle of a scaled
c     hessian (approximation) and a nonzero scaled gradient vector,
c     this subroutine computes a goldfeld-quandt-trotter step of
c     approximate length v(radius) by the more-hebden technique.  in
c     other words, step is computed to (approximately) minimize
c     psi(step) = (g**t)*step + 0.5*(step**t)*h*step  such that the
c     2-norm of d*step is at most (approximately) v(radius), where
c     g  is the gradient,  h  is the hessian, and  d  is a diagonal
c     scale matrix whose diagonal is stored in the parameter d.
c     (gqtstp assumes  dig = d**-1 * g  and  dihdi = d**-1 * h * d**-1.)
c     if g = 0, however, step = 0 is returned (even at a saddle point).
c
c  ***  parameter description  ***
c
c     d (in)  = the scale vector, i.e. the diagonal of the scale
c              matrix  d  mentioned above under purpose.
c   dig (in)  = the scaled gradient vector, d**-1 * g.  if g = 0, then
c              step = 0  and  v(stppar) = 0  are returned.
c dihdi (in)  = lower triangle of the scaled hessian (approximation),
c              i.e., d**-1 * h * d**-1, stored compactly by rows., i.e.,
c              in the order (1,1), (2,1), (2,2), (3,1), (3,2), etc.
c    ka (i/o) = the number of hebden iterations (so far) taken to deter-
c              mine step.  ka .lt. 0 on input means this is the first
c              attempt to determine step (for the present dig and dihdi)
c              -- ka is initialized to 0 in this case.  output with
c              ka = 0  (or v(stppar) = 0)  means  step = -(h**-1)*g.
c     l (i/o) = workspace of length p*(p+1)/2 for cholesky factors.
c     p (in)  = number of parameters -- the hessian is a  p x p  matrix.
c  step (i/o) = the step computed.
c     v (i/o) contains various constants and variables described below.
c     w (i/o) = workspace of length 4*p + 6.
c
c  ***  entries in v  ***
c
c v(dgnorm) (i/o) = 2-norm of (d**-1)*g.
c v(dstnrm) (output) = 2-norm of d*step.
c v(dst0)   (i/o) = 2-norm of d*(h**-1)*g (for pos. def. h only), or
c             overestimate of smallest eigenvalue of (d**-1)*h*(d**-1).
c v(epslon) (in)  = max. rel. error allowed for psi(step).  for the
c             step returned, psi(step) will exceed its optimal value
c             by less than -v(epslon)*psi(step).  suggested value = 0.1.
c v(gtstep) (out) = inner product between g and step.
c v(nreduc) (out) = psi(-(h**-1)*g) = psi(newton step)  (for pos. def.
c             h only -- v(nreduc) is set to zero otherwise).
c v(phmnfc) (in)  = tol. (together with v(phmxfc)) for accepting step
c             (more*s sigma).  the error v(dstnrm) - v(radius) must lie
c             between v(phmnfc)*v(radius) and v(phmxfc)*v(radius).
c v(phmxfc) (in)  (see v(phmnfc).)
c             suggested values -- v(phmnfc) = -0.25, v(phmxfc) = 0.5.
c v(preduc) (out) = psi(step) = predicted obj. func. reduction for step.
c v(radius) (in)  = radius of current (scaled) trust region.
c v(rad0)   (i/o) = value of v(radius) from previous call.
c v(stppar) (i/o) is normally the marquardt parameter, i.e. the alpha
c             described below under algorithm notes.  if h + alpha*d**2
c             (see algorithm notes) is (nearly) singular, however,
c             then v(stppar) = -alpha.
c
c  ***  usage notes  ***
c
c     if it is desired to recompute step using a different value of
c     v(radius), then this routine may be restarted by calling it
c     with all parameters unchanged except v(radius).  (this explains
c     why step and w are listed as i/o).  on an intiial call (one with
c     ka .lt. 0), step and w need not be initialized and only compo-
c     nents v(epslon), v(stppar), v(phmnfc), v(phmxfc), v(radius), and
c     v(rad0) of v must be initialized.  to compute step from a saddle
c     point (where the true gradient vanishes and h has a negative
c     eigenvalue), a nonzero g with small components should be passed.
c
c  ***  application and usage restrictions  ***
c
c     this routine is called as part of the nl2sol (nonlinear least-
c     squares) package (ref. 1), but it could be used in solving any
c     unconstrained minimization problem.
c
c  ***  algorithm notes  ***
c
c        the desired g-q-t step (ref. 2, 3, 4) satisfies
c     (h + alpha*d**2)*step = -g  for some nonnegative alpha such that
c     h + alpha*d**2 is positive semidefinite.  alpha and step are
c     computed by a scheme analogous to the one described in ref. 5.
c     estimates of the smallest and largest eigenvalues of the hessian
c     are obtained from the gerschgorin circle theorem enhanced by a
c     simple form of the scaling described in ref. 6.  cases in which
c     h + alpha*d**2 is nearly (or exactly) singular are handled by
c     the technique discussed in ref. 2.  in these cases, a step of
c     (exact) length v(radius) is returned for which psi(step) exceeds
c     its optimal value by less than -v(epslon)*psi(step).
c
c  ***  functions and subroutines called  ***
c
c dotprd - returns inner product of two vectors.
c litvmu - applies inverse-transpose of compact lower triang. matrix.
c livmul - applies inverse of compact lower triang. matrix.
c lsqrt  - finds cholesky factor (of compactly stored lower triang.).
c lsvmin - returns approx. to min. sing. value of lower triang. matrix.
c rmdcon - returns machine-dependent constants.
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
c 3.  goldfeld, s.m., quandt, r.e., and trotter, h.f. (1966),
c             maximization by quadratic hill-climbing, econometrica 34,
c             pp. 541-551.
c 4.  hebden, m.d. (1973), an algorithm for minimization using exact
c             second derivatives, report t.p. 515, theoretical physics
c             div., a.e.r.e. harwell, oxon., england.
c 5.  more, j.j. (1978), the levenberg-marquardt algorithm, implemen-
c             tation and theory, pp.105-116 of springer lecture notes
c             in mathematics no. 630, edited by g.a. watson, springer-
c             verlag, berlin and new york.
c 6.  varga, r.s. (1965), minimal gerschgorin sets, pacific j. math. 15,
c             pp. 719-729.
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
      logical restrt
      integer dggdmx, diag, diag0, dstsav, emax, emin, i, im1, inc, irc,
     1        j, k, kalim, k1, lk0, phipin, q, q0, uk0, x, x0
      double precision alphak, aki, akk, delta, dst, epso6, lk,
     1                 oldphi, phi, phimax, phimin, psifac, rad,
     2                 root, si, sk, sw, t, twopsi, t1, uk, wi
c
c     ***  constants  ***
      double precision dgxfac, epsfac, four, half, kappa, negone, one,
     1                 p001, six, three, two, zero
c
c  ***  intrinsic functions  ***
c/+
      double precision dabs, dmax1, dmin1, dsqrt
c/
c  ***  external functions and subroutines  ***
c
      external dotprd, litvmu, livmul, lsqrt, lsvmin, rmdcon, v2norm
      double precision dotprd, lsvmin, rmdcon, v2norm
c
c  ***  subscripts for v  ***
c
      integer dgnorm, dstnrm, dst0, epslon, gtstep, stppar, nreduc,
     1        phmnfc, phmxfc, preduc, radius, rad0
c/6
      data dgnorm/1/, dstnrm/2/, dst0/3/, epslon/19/,
     1     gtstep/4/, nreduc/6/, phmnfc/20/,
     2     phmxfc/21/, preduc/7/, radius/8/,
     3     rad0/9/, stppar/5/
c/7
c     parameter (dgnorm=1, dstnrm=2, dst0=3, epslon=19,
c    1     gtstep=4, nreduc=6, phmnfc=20,
c    2     phmxfc=21, preduc=7, radius=8,
c    3     rad0=9, stppar=5)
c/
c
c/6
      data epsfac/50.0d+0/, four/4.0d+0/, half/0.5d+0/,
     1     kappa/2.0d+0/, negone/-1.0d+0/, one/1.0d+0/, p001/1.0d-3/,
     2     six/6.0d+0/, three/3.0d+0/, two/2.0d+0/, zero/0.0d+0/
c/7
c     parameter (epsfac=50.0d+0, four=4.0d+0, half=0.5d+0,
c    1     kappa=2.0d+0, negone=-1.0d+0, one=1.0d+0, p001=1.0d-3,
c    2     six=6.0d+0, three=3.0d+0, two=2.0d+0, zero=0.0d+0)
c     save dgxfac
c/
      data dgxfac/0.d+0/
c
c  ***  body  ***
c
c     ***  store largest abs. entry in (d**-1)*h*(d**-1) at w(dggdmx).
      dggdmx = p + 1
c     ***  store gerschgorin over- and underestimates of the largest
c     ***  and smallest eigenvalues of (d**-1)*h*(d**-1) at w(emax)
c     ***  and w(emin) respectively.
      emax = dggdmx + 1
      emin = emax + 1
c     ***  for use in recomputing step, the final values of lk, uk, dst,
c     ***  and the inverse derivative of more*s phi at 0 (for pos. def.
c     ***  h) are stored in w(lk0), w(uk0), w(dstsav), and w(phipin)
c     ***  respectively.
      lk0 = emin + 1
      phipin = lk0 + 1
      uk0 = phipin + 1
      dstsav = uk0 + 1
c     ***  store diag of (d**-1)*h*(d**-1) in w(diag),...,w(diag0+p).
      diag0 = dstsav
      diag = diag0 + 1
c     ***  store -d*step in w(q),...,w(q0+p).
      q0 = diag0 + p
      q = q0 + 1
      rad = v(radius)
c     ***  phitol = max. error allowed in dst = v(dstnrm) = 2-norm of
c     ***  d*step.
      phimax = v(phmxfc) * rad
      phimin = v(phmnfc) * rad
c     ***  epso6 and psifac are used in checking for the special case
c     ***  of (nearly) singular h + alpha*d**2 (see ref. 2).
      psifac = two * v(epslon) / (three * (four * (v(phmnfc) + one) *
     1                       (kappa + one)  +  kappa  +  two) * rad**2)
c     ***  oldphi is used to detect limits of numerical accuracy.  if
c     ***  we recompute step and it does not change, then we accept it.
      oldphi = zero
      epso6 = v(epslon)/six
      irc = 0
      restrt = .false.
      kalim = ka + 50
c
c  ***  start or restart, depending on ka  ***
c
      if (ka .ge. 0) go to 310
c
c  ***  fresh start  ***
c
      k = 0
      uk = negone
      ka = 0
      kalim = 50
c
c     ***  store diag(dihdi) in w(diag0+1),...,w(diag0+p)  ***
c
      j = 0
      do 20 i = 1, p
         j = j + i
         k1 = diag0 + i
         w(k1) = dihdi(j)
 20      continue
c
c     ***  determine w(dggdmx), the largest element of dihdi  ***
c
      t1 = zero
      j = p * (p + 1) / 2
      do 30 i = 1, j
         t = dabs(dihdi(i))
         if (t1 .lt. t) t1 = t
 30      continue
      w(dggdmx) = t1
c
c  ***  try alpha = 0  ***
c
 40   call lsqrt(1, p, l, dihdi, irc)
      if (irc .eq. 0) go to 60
c        ***  indef. h -- underestimate smallest eigenvalue, use this
c        ***  estimate to initialize lower bound lk on alpha.
         j = irc*(irc+1)/2
         t = l(j)
         l(j) = one
         do 50 i = 1, irc
 50           w(i) = zero
         w(irc) = one
         call litvmu(irc, w, l, w)
         t1 = v2norm(irc, w)
         lk = -t / t1 / t1
         v(dst0) = -lk
         if (restrt) go to 210
         v(nreduc) = zero
         go to 70
c
c     ***  positive definite h -- compute unmodified newton step.  ***
 60   lk = zero
      call livmul(p, w(q), l, dig)
      v(nreduc) = half * dotprd(p, w(q), w(q))
      call litvmu(p, w(q), l, w(q))
      dst = v2norm(p, w(q))
      v(dst0) = dst
      phi = dst - rad
      if (phi .le. phimax) go to 280
      if (restrt) go to 210
c
c  ***  prepare to compute gerschgorin estimates of largest (and
c  ***  smallest) eigenvalues.  ***
c
 70   v(dgnorm) = v2norm(p, dig)
      if (v(dgnorm) .eq. zero) go to 450
      k = 0
      do 100 i = 1, p
         wi = zero
         if (i .eq. 1) go to 90
         im1 = i - 1
         do 80 j = 1, im1
              k = k + 1
              t = dabs(dihdi(k))
              wi = wi + t
              w(j) = w(j) + t
 80           continue
 90      w(i) = wi
         k = k + 1
 100     continue
c
c  ***  (under-)estimate smallest eigenvalue of (d**-1)*h*(d**-1)  ***
c
      k = 1
      t1 = w(diag) - w(1)
      if (p .le. 1) go to 120
      do 110 i = 2, p
         j = diag0 + i
         t = w(j) - w(i)
         if (t .ge. t1) go to 110
              t1 = t
              k = i
 110     continue
c
 120  sk = w(k)
      j = diag0 + k
      akk = w(j)
      k1 = k*(k-1)/2 + 1
      inc = 1
      t = zero
      do 150 i = 1, p
         if (i .eq. k) go to 130
         aki = dabs(dihdi(k1))
         si = w(i)
         j = diag0 + i
         t1 = half * (akk - w(j) + si - aki)
         t1 = t1 + dsqrt(t1*t1 + sk*aki)
         if (t .lt. t1) t = t1
         if (i .lt. k) go to 140
 130     inc = i
 140     k1 = k1 + inc
 150     continue
c
      w(emin) = akk - t
      uk = v(dgnorm)/rad - w(emin)
c
c  ***  compute gerschgorin (over-)estimate of largest eigenvalue  ***
c
      k = 1
      t1 = w(diag) + w(1)
      if (p .le. 1) go to 170
      do 160 i = 2, p
         j = diag0 + i
         t = w(j) + w(i)
         if (t .le. t1) go to 160
              t1 = t
              k = i
 160     continue
c
 170  sk = w(k)
      j = diag0 + k
      akk = w(j)
      k1 = k*(k-1)/2 + 1
      inc = 1
      t = zero
      do 200 i = 1, p
         if (i .eq. k) go to 180
         aki = dabs(dihdi(k1))
         si = w(i)
         j = diag0 + i
         t1 = half * (w(j) + si - aki - akk)
         t1 = t1 + dsqrt(t1*t1 + sk*aki)
         if (t .lt. t1) t = t1
         if (i .lt. k) go to 190
 180     inc = i
 190     k1 = k1 + inc
 200     continue
c
      w(emax) = akk + t
      lk = dmax1(lk, v(dgnorm)/rad - w(emax))
c
c     ***  alphak = current value of alpha (see alg. notes above).  we
c     ***  use more*s scheme for initializing it.
      alphak = dabs(v(stppar)) * v(rad0)/rad
c
      if (irc .ne. 0) go to 210
c
c  ***  compute l0 for positive definite h  ***
c
      call livmul(p, w, l, w(q))
      t = v2norm(p, w)
      w(phipin) = dst / t / t
      lk = dmax1(lk, phi*w(phipin))
c
c  ***  safeguard alphak and add alphak*i to (d**-1)*h*(d**-1)  ***
c
 210  ka = ka + 1
      if (-v(dst0) .ge. alphak .or. alphak .lt. lk .or. alphak .ge. uk)
     1                      alphak = uk * dmax1(p001, dsqrt(lk/uk))
      k = 0
      do 220 i = 1, p
         k = k + i
         j = diag0 + i
         dihdi(k) = w(j) + alphak
 220     continue
c
c  ***  try computing cholesky decomposition  ***
c
      call lsqrt(1, p, l, dihdi, irc)
      if (irc .eq. 0) go to 250
c
c  ***  (d**-1)*h*(d**-1) + alphak*i  is indefinite -- overestimate
c  ***  smallest eigenvalue for use in updating lk  ***
c
      j = (irc*(irc+1))/2
      t = l(j)
      l(j) = one
      do 230 i = 1, irc
 230     w(i) = zero
      w(irc) = one
      call litvmu(irc, w, l, w)
      t1 = v2norm(irc, w)
      lk = alphak - t/t1/t1
      v(dst0) = -lk
      go to 210
c
c  ***  alphak makes (d**-1)*h*(d**-1) positive definite.
c  ***  compute q = -d*step, check for convergence.  ***
c
 250  call livmul(p, w(q), l, dig)
      call litvmu(p, w(q), l, w(q))
      dst = v2norm(p, w(q))
      phi = dst - rad
      if (phi .le. phimax .and. phi .ge. phimin) go to 290
      if (phi .eq. oldphi) go to 290
      oldphi = phi
      if (phi .gt. zero) go to 260
c        ***  check for the special case of  h + alpha*d**2  (nearly)
c        ***  singular.  delta is .ge. the smallest eigenvalue of
c        ***  (d**-1)*h*(d**-1) + alphak*i.
         if (v(dst0) .gt. zero) go to 260
         delta = alphak + v(dst0)
         twopsi = alphak*dst*dst + dotprd(p, dig, w(q))
         if (delta .lt. psifac*twopsi) go to 270
c
c  ***  unacceptable alphak -- update lk, uk, alphak  ***
c
 260  if (ka .ge. kalim) go to 290
      call livmul(p, w, l, w(q))
      t1 = v2norm(p, w)
c     ***  the following dmin1 is necessary because of restarts  ***
      if (phi .lt. zero) uk = dmin1(uk, alphak)
      alphak = alphak  +  (phi/t1) * (dst/t1) * (dst/rad)
      lk = dmax1(lk, alphak)
      go to 210
c
c  ***  decide how to handle (nearly) singular h + alpha*d**2  ***
c
c     ***  if not yet available, obtain machine dependent value dgxfac.
 270  if (dgxfac .eq. zero) dgxfac = epsfac * rmdcon(3)
c
c     ***  now decide.  ***
      if (delta .gt. dgxfac*w(dggdmx)) go to 350
c        ***  delta is so small we cannot handle the special case in
c        ***  the available arithmetic.  accept step as it is.
         go to 290
c
c  ***  acceptable step on first try  ***
c
 280  alphak = zero
c
c  ***  successful step in general.  compute step = -(d**-1)*q  ***
c
 290  do 300 i = 1, p
         j = q0 + i
         step(i) = -w(j)/d(i)
 300     continue
      v(gtstep) = -dotprd(p, dig, w(q))
      v(preduc) = half * (dabs(alphak)*dst*dst - v(gtstep))
      go to 430
c
c
c  ***  restart with new radius  ***
c
 310  if (v(dst0) .le. zero .or. v(dst0) - rad .gt. phimax) go to 330
c
c     ***  prepare to return newton step  ***
c
         restrt = .true.
         ka = ka + 1
         k = 0
         do 320 i = 1, p
              k = k + i
              j = diag0 + i
              dihdi(k) = w(j)
 320          continue
         uk = negone
         go to 40
c
 330  if (ka .eq. 0) go to 60
c
      dst = w(dstsav)
      alphak = dabs(v(stppar))
      phi = dst - rad
      t = v(dgnorm)/rad
      if (rad .gt. v(rad0)) go to 340
c
c        ***  smaller radius  ***
         uk = t - w(emin)
         lk = zero
         if (alphak .gt. zero) lk = w(lk0)
         lk = dmax1(lk, t - w(emax))
         if (v(dst0) .gt. zero) lk = dmax1(lk, (v(dst0)-rad)*w(phipin))
         go to 260
c
c     ***  bigger radius  ***
 340  uk = t - w(emin)
      if (alphak .gt. zero) uk = dmin1(uk, w(uk0))
      lk = dmax1(zero, -v(dst0), t - w(emax))
      if (v(dst0) .gt. zero) lk = dmax1(lk, (v(dst0)-rad)*w(phipin))
      go to 260
c
c  ***  handle (nearly) singular h + alpha*d**2  ***
c
c     ***  negate alphak to indicate special case  ***
 350  alphak = -alphak
c     ***  allocate storage for scratch vector x  ***
      x0 = q0 + p
      x = x0 + 1
c
c  ***  use inverse power method with start from lsvmin to obtain
c  ***  approximate eigenvector corresponding to smallest eigenvalue
c  ***  of (d**-1)*h*(d**-1).
c
      delta = kappa*delta
      t = lsvmin(p, l, w(x), w)
c
      k = 0
c     ***  normalize w  ***
 360  do 370 i = 1, p
 370     w(i) = t*w(i)
c     ***  complete current inv. power iter. -- replace w by (l**-t)*w.
      call litvmu(p, w, l, w)
      t1 = one/v2norm(p, w)
      t = t1*t
      if (t .le. delta) go to 390
      if (k .gt. 30) go to 290
      k = k + 1
c     ***  start next inv. power iter. by storing normalized w in x.
      do 380 i = 1, p
         j = x0 + i
         w(j) = t1*w(i)
 380     continue
c     ***  compute w = (l**-1)*x.
      call livmul(p, w, l, w(x))
      t = one/v2norm(p, w)
      go to 360
c
 390  do 400 i = 1, p
 400     w(i) = t1*w(i)
c
c  ***  now w is the desired approximate (unit) eigenvector and
c  ***  t*x = ((d**-1)*h*(d**-1) + alphak*i)*w.
c
      sw = dotprd(p, w(q), w)
      t1 = (rad + dst) * (rad - dst)
      root = dsqrt(sw*sw + t1)
      if (sw .lt. zero) root = -root
      si = t1 / (sw + root)
c     ***  accept current step if adding si*w would lead to a
c     ***  further relative reduction in psi of less than v(epslon)/3.
      v(preduc) = half*twopsi
      t1 = zero
      t = si*(alphak*sw - half*si*(alphak + t*dotprd(p,w(x),w)))
      if (t .lt. epso6*twopsi) go to 410
         v(preduc) = v(preduc) + t
         dst = rad
         t1 = -si
 410  do 420 i = 1, p
         j = q0 + i
         w(j) = t1*w(i) - w(j)
         step(i) = w(j) / d(i)
 420     continue
      v(gtstep) = dotprd(p, dig, w(q))
c
c  ***  save values for use in a possible restart  ***
c
 430  v(dstnrm) = dst
      v(stppar) = alphak
      w(lk0) = lk
      w(uk0) = uk
      v(rad0) = rad
      w(dstsav) = dst
c
c     ***  restore diagonal of dihdi  ***
c
      j = 0
      do 440 i = 1, p
         j = j + i
         k = diag0 + i
         dihdi(j) = w(k)
 440     continue
      go to 999
c
c  ***  special case -- g = 0  ***
c
 450  v(stppar) = zero
      v(preduc) = zero
      v(dstnrm) = zero
      v(gtstep) = zero
      do 460 i = 1, p
 460     step(i) = zero
c
 999  return
c
c  ***  last card of gqtstp follows  ***
      end
      subroutine itsmry(d, iv, p, v, x)
c
c  ***  print nl2sol (version 2.2) iteration summary  ***
c
c  ***  parameter declarations  ***
c
      integer iv(1), p
      double precision d(p), v(1), x(p)
c     dimension iv(*), v(*)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  local variables  ***
c
      integer cov1, g1, i, ii, iv1, i1, j, m, nf, ng, ol, pu
c/6
      real model1(6), model2(6)
c/7
c     character*4 model1(6), model2(6)
c/
      double precision nreldf, oldf, preldf, reldf, zero
c
c  ***  intrinsic functions  ***
c/+
      integer iabs
c/
c  ***  no external functions or subroutines  ***
c
c  ***  subscripts for iv and v  ***
c
      integer covmat, covprt, covreq, dstnrm, f, fdif, f0, g,
     1        needhd, nfcall, nfcov, ngcov, ngcall, niter, nreduc,
     2        outlev, preduc, prntit, prunit, reldx, size, solprt,
     3        statpr, stppar, sused, x0prt
c
c  ***  iv subscript values  ***
c
c/6
      data covmat/26/, covprt/14/, g/28/, covreq/15/,
     1     needhd/39/, nfcall/6/, nfcov/40/, ngcov/41/,
     2     ngcall/30/, niter/31/, outlev/19/, prntit/48/,
     3     prunit/21/, solprt/22/, statpr/23/, sused/57/,
     4     x0prt/24/
c/7
c     parameter (covmat=26, covprt=14, g=28, covreq=15,
c    1     needhd=39, nfcall=6, nfcov=40, ngcov=41,
c    2     ngcall=30, niter=31, outlev=19, prntit=48,
c    3     prunit=21, solprt=22, statpr=23, sused=57,
c    4     x0prt=24)
c/
c
c  ***  v subscript values  ***
c
c/6
      data dstnrm/2/, f/10/, f0/13/, fdif/11/, nreduc/6/,
     1     preduc/7/, reldx/17/, size/47/, stppar/5/
c/7
c     parameter (dstnrm=2, f=10, f0=13, fdif=11, nreduc=6,
c    1     preduc=7, reldx=17, size=47, stppar=5)
c/
c
c/6
      data zero/0.d+0/
c/7
c     parameter (zero=0.d+0)
c/
c/6
      data model1(1)/4h    /, model1(2)/4h    /, model1(3)/4h    /,
     1     model1(4)/4h    /, model1(5)/4h  g /, model1(6)/4h  s /,
     2     model2(1)/4h g  /, model2(2)/4h s  /, model2(3)/4hg-s /,
     3     model2(4)/4hs-g /, model2(5)/4h-s-g/, model2(6)/4h-g-s/
c/7
c     data model1/'    ','    ','    ','    ','  g ','  s '/,
c    1     model2/' g  ',' s  ','g-s ','s-g ','-s-g','-g-s'/
c/
c
c-----------------------------------------------------------------------
c
      pu = iv(prunit)
      if (pu .eq. 0) go to 999
      iv1 = iv(1)
      ol = iv(outlev)
      if (iv1 .lt. 2 .or. iv1 .gt. 15) go to 140
      if (ol .eq. 0) go to 20
      if (iv1 .ge. 12) go to 20
      if (iv1 .ge. 10 .and. iv(prntit) .eq. 0) go to 20
      if (iv1 .gt. 2) go to 10
         iv(prntit) = iv(prntit) + 1
         if (iv(prntit) .lt. iabs(ol)) go to 999
 10   nf = iv(nfcall) - iabs(iv(nfcov))
      iv(prntit) = 0
      reldf = zero
      preldf = zero
      oldf = v(f0)
      if (oldf .le. zero) go to 12
         reldf = v(fdif) / oldf
         preldf = v(preduc) / oldf
 12   if (ol .gt. 0) go to 15
c
c        ***  print short summary line  ***
c
         if (iv(needhd) .eq. 1) write(pu, 1010)
 1010 format(12h0   it    nf,6x,1hf,8x,5hreldf,6x,6hpreldf,5x,5hreldx)
         iv(needhd) = 0
         write(pu,1017) iv(niter), nf, v(f), reldf, preldf, v(reldx)
         go to 20
c
c     ***  print long summary line  ***
c
 15   if (iv(needhd) .eq. 1) write(pu,1015)
 1015 format(12h0   it    nf,6x,1hf,8x,5hreldf,6x,6hpreldf,5x,5hreldx,
     1       4x,15hmodel    stppar,6x,4hsize,6x,6hd*step,5x,7hnpreldf)
      iv(needhd) = 0
      m = iv(sused)
      nreldf = zero
      if (oldf .gt. zero) nreldf = v(nreduc) / oldf
      write(pu,1017) iv(niter), nf, v(f), reldf, preldf, v(reldx),
     1               model1(m), model2(m), v(stppar), v(size),
     2               v(dstnrm), nreldf
 1017 format(1x,i5,i6,4d11.3,a3,a4,4d11.3)
c
 20   go to (999,999,30,35,40,45,50,60,70,80,90,150,110,120,130), iv1
c
 30   write(pu,1030)
 1030 format(26h0***** x-convergence *****)
      go to 180
c
 35   write(pu,1035)
 1035 format(42h0***** relative function convergence *****)
      go to 180
c
 40   write(pu,1040)
 1040 format(49h0***** x- and relative function convergence *****)
      go to 180
c
 45   write(pu,1045)
 1045 format(42h0***** absolute function convergence *****)
      go to 180
c
 50   write(pu,1050)
 1050 format(33h0***** singular convergence *****)
      go to 180
c
 60   write(pu,1060)
 1060 format(30h0***** false convergence *****)
      go to 180
c
 70   write(pu,1070)
 1070 format(38h0***** function evaluation limit *****)
      go to 180
c
 80   write(pu,1080)
 1080 format(28h0***** iteration limit *****)
      go to 180
c
 90   write(pu,1090)
 1090 format(18h0***** stopx *****)
      go to 180
c
 110  write(pu,1100)
 1100 format(45h0***** initial sum of squares overflows *****)
c
      go to 150
c
 120  write(pu,1120)
 1120 format(37h0***** bad parameters to assess *****)
      go to 999
c
 130  write(pu,1130)
 1130 format(36h0***** j could not be computed *****)
      if (iv(niter) .gt. 0) go to 190
      go to 150
c
 140  write(pu,1140) iv1
 1140 format(14h0***** iv(1) =,i5,6h *****)
      go to 999
c
c  ***  initial call on itsmry  ***
c
 150  if (iv(x0prt) .ne. 0) write(pu,1150) (i, x(i), d(i), i = 1, p)
 1150 format(23h0    i     initial x(i),7x,4hd(i)//(1x,i5,d17.6,d14.3))
      if (iv1 .ge. 13) go to 999
      iv(needhd) = 0
      iv(prntit) = 0
      if (ol .eq. 0) go to 999
      if (ol .lt. 0) write(pu,1010)
      if (ol .gt. 0) write(pu,1015)
      write(pu,1160) v(f)
 1160 format(12h0    0     1,d11.3,11x,d11.3)
      go to 999
c
c  ***  print various information requested on solution  ***
c
 180  iv(needhd) = 1
      if (iv(statpr) .eq. 0) go to 190
         oldf = v(f0)
         preldf = zero
         nreldf = zero
         if (oldf .le. zero) go to 185
              preldf = v(preduc) / oldf
              nreldf = v(nreduc) / oldf
 185     nf = iv(nfcall) - iv(nfcov)
         ng = iv(ngcall) - iv(ngcov)
         write(pu,1180) v(f), v(reldx), nf, ng, preldf, nreldf
 1180 format(9h0function,d17.6,8h   reldx,d20.6/12h func. evals,
     1   i8,9x,11hgrad. evals,i8/7h preldf,d19.6,3x,7hnpreldf,d18.6)
c
         if (iv(nfcov) .gt. 0) write(pu,1185) iv(nfcov)
 1185    format(1h0,i4,34h extra func. evals for covariance.)
         if (iv(ngcov) .gt. 0) write(pu,1186) iv(ngcov)
 1186    format(1x,i4,34h extra grad. evals for covariance.)
c
 190  if (iv(solprt) .eq. 0) go to 210
         iv(needhd) = 1
         g1 = iv(g)
         write(pu,1190)
 1190 format(22h0    i      final x(i),8x,4hd(i),10x,4hg(i)/)
         do 200 i = 1, p
              write(pu,1200) i, x(i), d(i), v(g1)
              g1 = g1 + 1
 200          continue
 1200    format(1x,i5,d17.6,2d14.3)
c
 210  if (iv(covprt) .eq. 0) go to 999
      cov1 = iv(covmat)
      iv(needhd) = 1
      if (cov1) 220, 230, 240
 220  if (-1 .eq. cov1) write(pu,1220)
 1220 format(43h0++++++ indefinite covariance matrix ++++++)
      if (-2 .eq. cov1) write(pu,1225)
 1225 format(52h0++++++ oversize steps in computing covariance +++++)
      go to 999
c
 230  write(pu,1230)
 1230 format(45h0++++++ covariance matrix not computed ++++++)
      go to 999
c
 240  i = iabs(iv(covreq))
      if (i .le. 1) write(pu,1241)
 1241 format(48h0covariance = scale * h**-1 * (j**t * j) * h**-1/)
      if (i .eq. 2) write(pu,1242)
 1242 format(27h0covariance = scale * h**-1/)
      if (i .ge. 3) write(pu,1243)
 1243 format(36h0covariance = scale * (j**t * j)**-1/)
      ii = cov1 - 1
      if (ol .le. 0) go to 260
      do 250 i = 1, p
         i1 = ii + 1
         ii = ii + i
         write(pu,1250) i, (v(j), j = i1, ii)
 250     continue
 1250 format(4h row,i3,2x,9d12.4/(9x,9d12.4))
      go to 999
c
 260  do 270 i = 1, p
         i1 = ii + 1
         ii = ii + i
         write(pu,1270) i, (v(j), j = i1, ii)
 270     continue
 1270 format(4h row,i3,2x,5d12.4/(9x,5d12.4))
c
 999  return
c  ***  last card of itsmry follows  ***
      end
      subroutine linvrt(n, lin, l)
c
c  ***  compute  lin = l**-1,  both  n x n  lower triang. stored   ***
c  ***  compactly by rows.  lin and l may share the same storage.  ***
c
c  ***  parameters  ***
c
      integer n
      double precision l(1), lin(1)
c     dimension l(n*(n+1)/2), lin(n*(n+1)/2)
c
c  ***  local variables  ***
c
      integer i, ii, im1, jj, j0, j1, k, k0, np1
      double precision one, t, zero
c/6
      data one/1.d+0/, zero/0.d+0/
c/7
c     parameter (one=1.d+0, zero=0.d+0)
c/
c
c  ***  body  ***
c
      np1 = n + 1
      j0 = n*(np1)/2
      do 30 ii = 1, n
         i = np1 - ii
         lin(j0) = one/l(j0)
         if (i .le. 1) go to 999
         j1 = j0
         im1 = i - 1
         do 20 jj = 1, im1
              t = zero
              j0 = j1
              k0 = j1 - jj
              do 10 k = 1, jj
                   t = t - l(k0)*lin(j0)
                   j0 = j0 - 1
                   k0 = k0 + k - i
 10                continue
              lin(j0) = t/l(k0)
 20           continue
         j0 = j0 - 1
 30      continue
 999  return
c  ***  last card of linvrt follows  ***
      end
      subroutine litvmu(n, x, l, y)
c
c  ***  solve  (l**t)*x = y,  where  l  is an  n x n  lower triangular
c  ***  matrix stored compactly by rows.  x and y may occupy the same
c  ***  storage.  ***
c
      integer n
      double precision x(n), l(1), y(n)
      integer i, ii, ij, im1, i0, j, np1
      double precision xi, zero
c/6
      data zero/0.d+0/
c/7
c     parameter (zero=0.d+0)
c/
c
      do 10 i = 1, n
 10      x(i) = y(i)
      np1 = n + 1
      i0 = n*(n+1)/2
      do 30 ii = 1, n
         i = np1 - ii
         xi = x(i)/l(i0)
         x(i) = xi
         if (i .le. 1) go to 999
         i0 = i0 - i
         if (xi .eq. zero) go to 30
         im1 = i - 1
         do 20 j = 1, im1
              ij = i0 + j
              x(j) = x(j) - xi*l(ij)
 20           continue
 30      continue
 999  return
c  ***  last card of litvmu follows  ***
      end
      subroutine livmul(n, x, l, y)
c
c  ***  solve  l*x = y, where  l  is an  n x n  lower triangular
c  ***  matrix stored compactly by rows.  x and y may occupy the same
c  ***  storage.  ***
c
      integer n
      double precision x(n), l(1), y(n)
      external dotprd
      double precision dotprd
      integer i, j, k
      double precision t, zero
c/6
      data zero/0.d+0/
c/7
c     parameter (zero=0.d+0)
c/
c
      do 10 k = 1, n
         if (y(k) .ne. zero) go to 20
         x(k) = zero
 10      continue
      go to 999
 20   j = k*(k+1)/2
      x(k) = y(k) / l(j)
      if (k .ge. n) go to 999
      k = k + 1
      do 30 i = k, n
         t = dotprd(i-1, l(j+1), x)
         j = j + i
         x(i) = (y(i) - t)/l(j)
 30      continue
 999  return
c  ***  last card of livmul follows  ***
      end
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
      data dgnorm/1/, dstnrm/2/, dst0/3/, epslon/19/,
     1     gtstep/4/, nreduc/6/, phmnfc/20/,
     2     phmxfc/21/, preduc/7/, radius/8/,
     3     rad0/9/, stppar/5/
c/7
c     parameter (dgnorm=1, dstnrm=2, dst0=3, epslon=19,
c    1     gtstep=4, nreduc=6, phmnfc=20,
c    2     phmxfc=21, preduc=7, radius=8,
c    3     rad0=9, stppar=5)
c/
c
c/6
      data dfac/256.d+0/, eight/8.d+0/, half/0.5d+0/, negone/-1.d+0/,
     1     one/1.d+0/, p001/1.d-3/, three/3.d+0/, ttol/2.5d+0/,
     2     zero/0.d+0/
c/7
c     parameter (dfac=256.d+0, eight=8.d+0, half=0.5d+0, negone=-1.d+0,
c    1     one=1.d+0, p001=1.d-3, three=3.d+0, ttol=2.5d+0,
c    2     zero=0.d+0)
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
      subroutine lsqrt(n1, n, l, a, irc)
c
c  ***  compute rows n1 through n of the cholesky factor  l  of
c  ***  a = l*(l**t),  where  l  and the lower triangle of  a  are both
c  ***  stored compactly by rows (and may occupy the same storage).
c  ***  irc = 0 means all went well.  irc = j means the leading
c  ***  principal  j x j  submatrix of  a  is not positive definite --
c  ***  and  l(j*(j+1)/2)  contains the (nonpos.) reduced j-th diagonal.
c
c  ***  parameters  ***
c
      integer n1, n, irc
      double precision l(1), a(1)
c     dimension l(n*(n+1)/2), a(n*(n+1)/2)
c
c  ***  local variables  ***
c
      integer i, ij, ik, im1, i0, j, jk, jm1, j0, k
      double precision t, td, zero
c
c  ***  intrinsic functions  ***
c/+
      double precision dsqrt
c/
c/6
      data zero/0.d+0/
c/7
c     parameter (zero=0.d+0)
c/
c
c  ***  body  ***
c
      i0 = n1 * (n1 - 1) / 2
      do 50 i = n1, n
         td = zero
         if (i .eq. 1) go to 40
         j0 = 0
         im1 = i - 1
         do 30 j = 1, im1
              t = zero
              if (j .eq. 1) go to 20
              jm1 = j - 1
              do 10 k = 1, jm1
                   ik = i0 + k
                   jk = j0 + k
                   t = t + l(ik)*l(jk)
 10                continue
 20           ij = i0 + j
              j0 = j0 + j
              t = (a(ij) - t) / l(j0)
              l(ij) = t
              td = td + t*t
 30           continue
 40      i0 = i0 + i
         t = a(i0) - td
         if (t .le. zero) go to 60
         l(i0) = dsqrt(t)
 50      continue
c
      irc = 0
      go to 999
c
 60   l(i0) = t
      irc = i
c
 999  return
c
c  ***  last card of lsqrt  ***
      end
      double precision function lsvmin(p, l, x, y)
c
c  ***  estimate smallest sing. value of packed lower triang. matrix l
c
c  ***  parameter declarations  ***
c
      integer p
      double precision l(1), x(p), y(p)
c     dimension l(p*(p+1)/2)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  purpose  ***
c
c     this function returns a good over-estimate of the smallest
c     singular value of the packed lower triangular matrix l.
c
c  ***  parameter description  ***
c
c  p (in)  = the order of l.  l is a  p x p  lower triangular matrix.
c  l (in)  = array holding the elements of  l  in row order, i.e.
c             l(1,1), l(2,1), l(2,2), l(3,1), l(3,2), l(3,3), etc.
c  x (out) if lsvmin returns a positive value, then x is a normalized
c             approximate left singular vector corresponding to the
c             smallest singular value.  this approximation may be very
c             crude.  if lsvmin returns zero, then some components of x
c             are zero and the rest retain their input values.
c  y (out) if lsvmin returns a positive value, then y = (l**-1)*x is an
c             unnormalized approximate right singular vector correspond-
c             ing to the smallest singular value.  this approximation
c             may be crude.  if lsvmin returns zero, then y retains its
c             input value.  the caller may pass the same vector for x
c             and y (nonstandard fortran usage), in which case y over-
c             writes x (for nonzero lsvmin returns).
c
c  ***  application and usage restrictions  ***
c
c     there are no usage restrictions.
c
c  ***  algorithm notes  ***
c
c     the algorithm is based on (1), with the additional provision that
c     lsvmin = 0 is returned if the smallest diagonal element of l
c     (in magnitude) is not more than the unit roundoff times the
c     largest.  the algorithm uses a random number generator proposed
c     in (4), which passes the spectral test with flying colors -- see
c     (2) and (3).
c
c  ***  subroutines and functions called  ***
c
c        v2norm - function, returns the 2-norm of a vector.
c
c  ***  references  ***
c
c     (1) cline, a., moler, c., stewart, g., and wilkinson, j.h.(1977),
c         an estimate for the condition number of a matrix, report
c         tm-310, applied math. div., argonne national laboratory.
c
c     (2) hoaglin, d.c. (1976), theoretical properties of congruential
c         random-number generators --  an empirical view,
c         memorandum ns-340, dept. of statistics, harvard univ.
c
c     (3) knuth, d.e. (1969), the art of computer programming, vol. 2
c         (seminumerical algorithms), addison-wesley, reading, mass.
c
c     (4) smith, c.s. (1971), multiplicative pseudo-random number
c         generators with prime modulus, j. assoc. comput. mach. 18,
c         pp. 586-593.
c
c  ***  history  ***
c
c     designed and coded by david m. gay (winter 1977/summer 1978).
c
c  ***  general  ***
c
c     this subroutine was written in connection with research
c     supported by the national science foundation under grants
c     mcs-7600324, dcr75-10143, 76-14311dss, and mcs76-11989.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  local variables  ***
c
      integer i, ii, ix, j, ji, jj, jjj, jm1, j0, pplus1
      double precision b, psj, sminus, splus, t, xminus, xplus
c
c  ***  constants  ***
c
      double precision half, one, r9973, zero
c
c  ***  intrinsic functions  ***
c/+
      integer mod
      real float
      double precision dabs
c/
c  ***  external functions and subroutines  ***
c
      external v2norm
      double precision v2norm
c
c/6
      data half/0.5d+0/, one/1.d+0/, r9973/9973.d+0/, zero/0.d+0/
c/7
c     parameter (half=0.5d+0, one=1.d+0, r9973=9973.d+0, zero=0.d+0)
c     save ix
c/
      data ix/2/
c
c  ***  body  ***
c
c  ***  first check whether to return lsvmin = 0 and initialize x  ***
c
      ii = 0
      do 10 i = 1, p
         x(i) = zero
         ii = ii + i
         if (l(ii) .eq. zero) go to 300
 10      continue
      if (mod(ix, 9973) .eq. 0) ix = 2
      pplus1 = p + 1
c
c  ***  solve (l**t)*x = b, where the components of b have randomly
c  ***  chosen magnitudes in (.5,1) with signs chosen to make x large.
c
c     do j = p to 1 by -1...
      do 100 jjj = 1, p
         j = pplus1 - jjj
c       ***  determine x(j) in this iteration. note for i = 1,2,...,j
c       ***  that x(i) holds the current partial sum for row i.
         ix = mod(3432*ix, 9973)
         b = half*(one + float(ix)/r9973)
         xplus = (b - x(j))
         xminus = (-b - x(j))
         splus = dabs(xplus)
         sminus = dabs(xminus)
         jm1 = j - 1
         j0 = j*jm1/2
         jj = j0 + j
         xplus = xplus/l(jj)
         xminus = xminus/l(jj)
         if (jm1 .eq. 0) go to 30
         do 20 i = 1, jm1
              ji = j0 + i
              splus = splus + dabs(x(i) + l(ji)*xplus)
              sminus = sminus + dabs(x(i) + l(ji)*xminus)
 20           continue
 30      if (sminus .gt. splus) xplus = xminus
         x(j) = xplus
c       ***  update partial sums  ***
         if (jm1 .eq. 0) go to 100
         do 40 i = 1, jm1
              ji = j0 + i
              x(i) = x(i) + l(ji)*xplus
 40           continue
 100     continue
c
c  ***  normalize x  ***
c
      t = one/v2norm(p, x)
      do 110 i = 1, p
 110     x(i) = t*x(i)
c
c  ***  solve l*y = x and return svmin = 1/twonorm(y)  ***
c
      do 200 j = 1, p
         psj = zero
         jm1 = j - 1
         j0 = j*jm1/2
         if (jm1 .eq. 0) go to 130
         do 120 i = 1, jm1
              ji = j0 + i
              psj = psj + l(ji)*y(i)
 120          continue
 130     jj = j0 + j
         y(j) = (x(j) - psj)/l(jj)
 200     continue
c
      lsvmin = one/v2norm(p, y)
      go to 999
c
 300  lsvmin = zero
 999  return
c  ***  last card of lsvmin follows  ***
      end
      subroutine ltsqar(n, a, l)
c
c  ***  set a to lower triangle of (l**t) * l  ***
c
c  ***  l = n x n lower triang. matrix stored rowwise.  ***
c  ***  a is also stored rowwise and may share storage with l.  ***
c
      integer n
      double precision a(1), l(1)
c     dimension a(n*(n+1)/2), l(n*(n+1)/2)
c
      integer i, ii, iim1, i1, j, k, m
      double precision lii, lj
c
      ii = 0
      do 50 i = 1, n
         i1 = ii + 1
         ii = ii + i
         m = 1
         if (i .eq. 1) go to 30
         iim1 = ii - 1
         do 20 j = i1, iim1
              lj = l(j)
              do 10 k = i1, j
                   a(m) = a(m) + lj*l(k)
                   m = m + 1
 10                continue
 20           continue
 30      lii = l(ii)
         do 40 j = i1, ii
 40           a(j) = lii * l(j)
 50      continue
c
 999  return
c  ***  last card of ltsqar follows  ***
      end
      subroutine parchk(iv, n, nn, p, v)
c
c  ***  check nl2sol (version 2.2) parameters, print changed values  ***
c
      integer iv(1), n, nn, p
      double precision v(1)
c     dimension iv(*), v(*)
c
      external dfault, rmdcon, vcopy
      double precision rmdcon
c dfault -- supplies dfault parameter values.
c rmdcon -- returns machine-dependent constants.
c vcopy  -- copies one vector to another.
c
c  ***  local variables  ***
c
      integer i, iv1, jtolp, k, l, m, nvdflt, pu
c/6
      real cngd(3), dflt(3), vn(2,27), which(3)
c/7
c     character*4 cngd(3), dflt(3), vn(2,27), which(3)
c/
      double precision big, machep, tiny, vk, vm(27), vx(27), zero
c
c  ***  iv and v subscripts  ***
c
      integer dtype, dtype0, d0init, epslon, inits, jtinit, jtol0,
     1        jtol1, oldn, oldnn, oldp, parprt, parsv1, prunit
c
c/6
      data nvdflt/27/, zero/0.d+0/
c/7
c     parameter (nvdflt=27, zero=0.d+0)
c/
c
c/6
      data dtype/16/, dtype0/29/, d0init/37/, epslon/19/,
     1     inits/25/, jtinit/39/, jtol0/86/, jtol1/87/,
     2     oldn/45/, oldnn/46/, oldp/47/, parprt/20/,
     3     parsv1/51/, prunit/21/
c/7
c     parameter (dtype=16, dtype0=29, d0init=37, epslon=19,
c    1     inits=25, jtinit=39, jtol0=86, jtol1=87,
c    2     oldn=45, oldnn=46, oldp=47, parprt=20,
c    3     parsv1=51, prunit=21)
c     save big, tiny
c/
c
      data big/0.d+0/, tiny/1.d+0/
c/6
      data vn(1,1),vn(2,1)/4hepsl,4hon../
      data vn(1,2),vn(2,2)/4hphmn,4hfc../
      data vn(1,3),vn(2,3)/4hphmx,4hfc../
      data vn(1,4),vn(2,4)/4hdecf,4hac../
      data vn(1,5),vn(2,5)/4hincf,4hac../
      data vn(1,6),vn(2,6)/4hrdfc,4hmn../
      data vn(1,7),vn(2,7)/4hrdfc,4hmx../
      data vn(1,8),vn(2,8)/4htune,4hr1../
      data vn(1,9),vn(2,9)/4htune,4hr2../
      data vn(1,10),vn(2,10)/4htune,4hr3../
      data vn(1,11),vn(2,11)/4htune,4hr4../
      data vn(1,12),vn(2,12)/4htune,4hr5../
      data vn(1,13),vn(2,13)/4hafct,4hol../
      data vn(1,14),vn(2,14)/4hrfct,4hol../
      data vn(1,15),vn(2,15)/4hxcto,4hl.../
      data vn(1,16),vn(2,16)/4hxfto,4hl.../
      data vn(1,17),vn(2,17)/4hlmax,4h0.../
      data vn(1,18),vn(2,18)/4hdltf,4hdj../
      data vn(1,19),vn(2,19)/4hd0in,4hit../
      data vn(1,20),vn(2,20)/4hdini,4ht.../
      data vn(1,21),vn(2,21)/4hjtin,4hit../
      data vn(1,22),vn(2,22)/4hdltf,4hdc../
      data vn(1,23),vn(2,23)/4hdfac,4h..../
      data vn(1,24),vn(2,24)/4hrlim,4hit../
      data vn(1,25),vn(2,25)/4hcosm,4hin../
      data vn(1,26),vn(2,26)/4hdelt,4ha0../
      data vn(1,27),vn(2,27)/4hfuzz,4h..../
c/7
c     data vn(1,1),vn(2,1)/'epsl','on..'/
c     data vn(1,2),vn(2,2)/'phmn','fc..'/
c     data vn(1,3),vn(2,3)/'phmx','fc..'/
c     data vn(1,4),vn(2,4)/'decf','ac..'/
c     data vn(1,5),vn(2,5)/'incf','ac..'/
c     data vn(1,6),vn(2,6)/'rdfc','mn..'/
c     data vn(1,7),vn(2,7)/'rdfc','mx..'/
c     data vn(1,8),vn(2,8)/'tune','r1..'/
c     data vn(1,9),vn(2,9)/'tune','r2..'/
c     data vn(1,10),vn(2,10)/'tune','r3..'/
c     data vn(1,11),vn(2,11)/'tune','r4..'/
c     data vn(1,12),vn(2,12)/'tune','r5..'/
c     data vn(1,13),vn(2,13)/'afct','ol..'/
c     data vn(1,14),vn(2,14)/'rfct','ol..'/
c     data vn(1,15),vn(2,15)/'xcto','l...'/
c     data vn(1,16),vn(2,16)/'xfto','l...'/
c     data vn(1,17),vn(2,17)/'lmax','0...'/
c     data vn(1,18),vn(2,18)/'dltf','dj..'/
c     data vn(1,19),vn(2,19)/'d0in','it..'/
c     data vn(1,20),vn(2,20)/'dini','t...'/
c     data vn(1,21),vn(2,21)/'jtin','it..'/
c     data vn(1,22),vn(2,22)/'dltf','dc..'/
c     data vn(1,23),vn(2,23)/'dfac','....'/
c     data vn(1,24),vn(2,24)/'rlim','it..'/
c     data vn(1,25),vn(2,25)/'cosm','in..'/
c     data vn(1,26),vn(2,26)/'delt','a0..'/
c     data vn(1,27),vn(2,27)/'fuzz','....'/
c/
c
      data vm(1)/1.0d-3/, vm(2)/-0.99d+0/, vm(3)/1.0d-3/, vm(4)/1.0d-2/,
     1     vm(5)/1.2d+0/, vm(6)/1.d-2/, vm(7)/1.2d+0/, vm(8)/0.d+0/,
     2     vm(9)/0.d+0/, vm(10)/1.d-3/, vm(11)/-1.d+0/, vm(15)/0.d+0/,
     3     vm(16)/0.d+0/, vm(19)/0.d+0/, vm(20)/-10.d+0/, vm(21)/0.d+0/,
     4     vm(23)/0.d+0/, vm(24)/1.d+10/, vm(27)/1.01d+0/
      data vx(1)/0.9d+0/, vx(2)/-1.d-3/, vx(3)/1.d+1/, vx(4)/0.8d+0/,
     1     vx(5)/1.d+2/, vx(6)/0.8d+0/, vx(7)/1.d+2/, vx(8)/0.5d+0/,
     2     vx(9)/0.5d+0/, vx(10)/1.d+0/, vx(11)/1.d+0/, vx(14)/0.1d+0/,
     3     vx(15)/1.d+0/, vx(16)/1.d+0/, vx(18)/1.d+0/, vx(22)/1.d+0/,
     4     vx(23)/1.d+0/, vx(25)/1.d+0/, vx(26)/1.d+0/, vx(27)/1.d+2/
c
c/6
      data cngd(1),cngd(2),cngd(3)/4h---c,4hhang,4hed v/,
     1     dflt(1),dflt(2),dflt(3)/4hnond,4hefau,4hlt v/
c/7
c     data cngd(1),cngd(2),cngd(3)/'---c','hang','ed v'/,
c    1     dflt(1),dflt(2),dflt(3)/'nond','efau','lt v'/
c/
c
c.......................................................................
c
      if (iv(1) .eq. 0) call dfault(iv, v)
      pu = iv(prunit)
      iv1 = iv(1)
      if (iv1 .ne. 12) go to 30
         if (nn .ge. n .and. n .ge. p .and. p .ge. 1) go to 20
              iv(1) = 16
              if (pu .ne. 0) write(pu,10) nn, n, p
 10           format(30h0///// bad nn, n, or p... nn =,i5,5h, n =,i5,
     1               5h, p =,i5)
              go to 999
 20      k = iv(21)
         call dfault(iv(21), v(33))
         iv(21) = k
         iv(dtype0) = iv(dtype+20)
         iv(oldn) = n
         iv(oldnn) = nn
         iv(oldp) = p
         which(1) = dflt(1)
         which(2) = dflt(2)
         which(3) = dflt(3)
         go to 80
 30   if (n .eq. iv(oldn) .and. nn .eq. iv(oldnn) .and. p .eq. iv(oldp))
     1                       go to 50
         iv(1) = 17
         if (pu .ne. 0) write(pu,40) iv(oldnn), iv(oldn), iv(oldp), nn,
     1                               n, p
 40      format(30h0///// (nn,n,p) changed from (,i5,1h,,i5,1h,,i3,
     1          6h) to (,i5,1h,,i5,1h,,i3,2h).)
         go to 999
c
 50   if (iv1 .le. 11 .and. iv1 .ge. 1) go to 70
         iv(1) = 50
         if (pu .ne. 0) write(pu,60) iv1
 60      format(15h0/////  iv(1) =,i5,28h should be between 0 and 12.)
         go to 999
c
 70   which(1) = cngd(1)
      which(2) = cngd(2)
      which(3) = cngd(3)
c
 80   if (big .gt. tiny) go to 90
         tiny = rmdcon(1)
         machep = rmdcon(3)
         big = rmdcon(6)
         vm(12) = machep
         vx(12) = big
         vm(13) = tiny
         vx(13) = big
         vm(14) = machep
         vm(17) = tiny
         vx(17) = big
         vm(18) = machep
         vx(19) = big
         vx(20) = big
         vx(21) = big
         vm(22) = machep
         vx(24) = rmdcon(5)
         vm(25) = machep
         vm(26) = machep
 90   m = 0
      if (iv(inits) .ge. 0 .and. iv(inits) .le. 2) go to 110
         m = 18
         if (pu .ne. 0) write(pu,100) iv(inits)
 100     format(25h0/////  inits... iv(25) =,i4,20h should be between 0,
     1          7h and 2.)
 110  k = epslon
      do 140 i = 1, nvdflt
         vk = v(k)
         if (vk .ge. vm(i) .and. vk .le. vx(i)) go to 130
              m = k
              if (pu .ne. 0) write(pu,120) vn(1,i), vn(2,i), k, vk,
     1                                    vm(i), vx(i)
 120          format(8h0/////  ,2a4,5h.. v(,i2,3h) =,d11.3,7h should,
     1               11h be between,d11.3,4h and,d11.3)
 130     k = k + 1
 140     continue
c
      if (iv1 .eq. 12 .and. v(jtinit) .gt. zero) go to 170
c
c  ***  check jtol values  ***
c
      jtolp = jtol0 + p
      do 160 i = jtol1, jtolp
         if (v(i) .gt. zero) go to 160
         k = i - jtol0
         if (pu .ne. 0) write(pu,150) k, i, v(i)
 150     format(12h0///// jtol(,i3,6h) = v(,i3,3h) =,d11.3,
     1          20h should be positive.)
         m = i
 160     continue
c
 170  if (m .eq. 0) go to 180
         iv(1) = m
         go to 999
c
 180  if (pu .eq. 0 .or. iv(parprt) .eq. 0) go to 999
      if (iv1 .ne. 12 .or. iv(inits) .eq. 0) go to 200
         m = 1
         write(pu,190) iv(inits)
 190     format(22h0nondefault values..../20h inits..... iv(25) =,i3)
 200  if (iv(dtype) .eq. iv(dtype0)) go to  210
         if (m .eq. 0) write(pu,215) which
         m = 1
         write(pu,205) iv(dtype)
 205     format(20h dtype..... iv(16) =,i3)
 210  k = epslon
      l = parsv1
      do 240 i = 1, nvdflt
         if (v(k) .eq. v(l)) go to 230
              if (m .eq. 0) write(pu,215) which
 215          format(1h0,3a4,9halues..../)
              m = 1
              write(pu,220) vn(1,i), vn(2,i), k, v(k)
 220          format(1x,2a4,5h.. v(,i2,3h) =,d15.7)
 230     k = k + 1
         l = l + 1
 240     continue
      iv(dtype0) = iv(dtype)
      call vcopy(nvdflt, v(parsv1), v(epslon))
      if (iv1 .ne. 12) go to 999
         if (v(jtinit) .gt. zero) go to 260
              jtolp = jtol0 + p
              write(pu,250) (v(i), i = jtol1, jtolp)
 250          format(24h0(initial) jtol array.../(1x,6d12.3))
 260     if (v(d0init) .gt. zero) go to 999
              k = jtol1 + p
              l = k + p - 1
              write(pu,270) (v(i), i = k, l)
 270          format(22h0(initial) d0 array.../1x,6d12.3)
c
 999  return
c  ***  last card of parchk follows  ***
      end
      subroutine qapply(nn, n, p, j, r, ierr)
c     *****parameters.
      integer nn, n, p, ierr
      double precision j(nn,p), r(n)
c
c     ..................................................................
c     ..................................................................
c
c     *****purpose.
c     this subroutine applies to r the orthogonal transformations
c     stored in j by qrfact
c
c     *****parameter description.
c     on input.
c
c        nn is the row dimension of the matrix j as declared in
c             the calling program dimension statement
c
c        n is the number of rows of j and the size of the vector r
c
c        p is the number of columns of j and the size of sigma
c
c        j contains on and below its diagonal the column vectors
c             u which determine the householder transformations
c             ident - u*u.transpose
c
c        r is the right hand side vector to which the orthogonal
c             transformations will be applied
c
c        ierr if non-zero indicates that not all the transformations
c             were successfully determined and only the first
c             abs(ierr) - 1 transformations will be used
c
c     on output.
c
c        r has been overwritten by its transformed image
c
c     *****application and usage restrictions.
c     none
c
c     *****algorithm notes.
c     the vectors u which determine the householder transformations
c     are normalized so that their 2-norm squared is 2.  the use of
c     these transformations here is in the spirit of (1).
c
c     *****subroutines and functions called.
c
c     dotprd - function, returns the inner product of vectors
c
c     *****references.
c     (1) businger, p. a., and golub, g. h. (1965), linear least squares
c        solutions by householder transformations, numer. math. 7,
c        pp. 269-276.
c
c     *****history.
c     designed by david m. gay, coded by stephen c. peters (winter 1977)
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
c     *****local variables.
      integer i, k, l, nl1
      double precision t
c     *****intrinsic functions.
c/+
      integer iabs
c/
c     *****functions.
      external dotprd
      double precision dotprd
c
      k = p
      if (ierr .ne. 0) k = iabs(ierr) - 1
      if ( k .eq. 0) go to 999
c
      do 20 l = 1, k
         nl1 = n - l + 1
         t = -dotprd(nl1, j(l,l), r(l))
c
         do 10 i = l, n
 10           r(i) = r(i) + t*j(i,l)
 20   continue
 999  return
c     .... last card of qapply .........................................
      end
      subroutine qrfact(nm,m,n,qr,alpha,ipivot,ierr,nopivk,sum)
c
c  ***  compute the qr decomposition of the matrix stored in qr  ***
c
c     *****parameters.
      integer nm,m,n,ipivot(n),ierr,nopivk
      double precision  qr(nm,n),alpha(n),sum(n)
c     *****local variables.
      integer i,j,jbar,k,k1,minum,mk1
      double precision  alphak,beta,qrkk,qrkmax,sigma,temp,ufeta,rktol,
     1        rktol1,sumj
c     *****functions.
c/+
      integer min0
      double precision  dabs,dsqrt
c/
      external dotprd, rmdcon, vaxpy, vscopy, v2norm
      double precision dotprd, rmdcon, v2norm
c dotprd... returns inner product of two vectors.
c rmdcon... returns machine-dependent constants.
c vaxpy... computes scalar times one vector plus another.
c vscopy... sets all elements of a vector to a scalar.
c v2norm... returns the 2-norm of a vector.
c
c     *****constants.
      double precision one, p01, p99, zero
c/6
      data one/1.0d+0/, p01/0.01d+0/, p99/0.99d+0/, zero/0.0d+0/
c/7
c     parameter (one=1.0d+0, p01=0.01d+0, p99=0.99d+0, zero=0.0d+0)
c     save rktol, ufeta
c/
c
c
c     ..................................................................
c     ..................................................................
c
c
c     *****purpose.
c
c     this subroutine does a qr-decomposition on the m x n matrix qr,
c        with an optionally modified column pivoting, and returns the
c        upper triangular r-matrix, as well as the orthogonal vectors
c        used in the transformations.
c
c     *****parameter description.
c     on input.
c
c        nm must be set to the row dimension of the two dimensional
c             array parameters as declared in the calling program
c             dimension statement.
c
c        m must be set to the number of rows in the matrix.
c
c        n must be set to the number of columns in the matrix.
c
c        qr contains the real rectangular matrix to be decomposed.
c
c     nopivk is used to control pivotting.  columns 1 through
c        nopivk will remain fixed in position.
c
c        sum is used for temporary storage for the subroutine.
c
c     on output.
c
c        qr contains the non-diagonal elements of the r-matrix
c             in the strict upper triangle. the vectors u, which
c             define the householder transformations   i - u*u-transp,
c             are in the columns of the lower triangle. these vectors u
c             are scaled so that the square of their 2-norm is 2.0.
c
c        alpha contains the diagonal elements of the r-matrix.
c
c        ipivot reflects the column pivoting performed on the input
c             matrix to accomplish the decomposition. the j-th
c             element of ipivot gives the column of the original
c             matrix which was pivoted into column j during the
c             decomposition.
c
c        ierr is set to.
c             0 for normal return,
c             k if no non-zero pivot could be found for the k-th
c                  transformation, or
c             -k for an error exit on the k-th thansformation.
c             if an error exit was taken, the first (k - 1)
c             transformations are correct.
c
c
c     *****applications and usage restrictions.
c     this may be used when solving linear least-squares problems --
c     see subroutine qr1 of rosepack.  it is called for this purpose
c     by llsqst in the nl2sol (nonlinear least-squares) package.
c
c     *****algorithm notes.
c     this version of qrfact tries to eliminate the occurrence of
c     underflows during the accumulation of inner products.  rktol1
c     is chosen below so as to insure that discarded terms have no
c     effect on the computed two-norms.
c
c     adapted from the algol routine solve (1).
c
c     *****references.
c     (1)     businger,p. and golub,g.h., linear least squares
c     solutions by housholder transformations, in wilkinson,j.h.
c     and reinsch,c.(eds.), handbook for automatic computation,
c     volume ii. linear algebra, springer-verlag, 111-118 (1971).
c     prepublished in numer.math. 7, 269-276 (1965).
c
c     *****history.
c     this amounts to the subroutine qr1 of rosepack with rktol1 used
c     in place of rktol below, with v2norm used to initialize (and
c     sometimes update) the sum array, and with calls on dotprd and
c     vaxpy in place of some loops.
c
c     *****general.
c
c     development of this program supported in part by
c     national science foundation grant gj-1154x3 and
c     national science foundation grant dcr75-08802
c     to national bureau of economic research, inc.
c
c
c
c     ..................................................................
c     ..................................................................
c
c
c     ..........  ufeta is the smallest positive floating point number
c        s.t. ufeta and -ufeta can both be represented.
c
c     ..........  rktol is the square root of the relative precision
c        of floating point arithmetic (machep).
      data rktol/0.d+0/, ufeta/0.d+0/
c     *****body of program.
      if (ufeta .gt. zero) go to 10
         ufeta = rmdcon(1)
         rktol = rmdcon(4)
   10 ierr = 0
      rktol1 = p01 * rktol
c
      do 20 j=1,n
         sum(j) = v2norm(m, qr(1,j))
         ipivot(j) = j
   20 continue
c
      minum = min0(m,n)
c
      do 120 k=1,minum
         mk1 = m - k + 1
c        ..........k-th householder transformation..........
         sigma = zero
         jbar = 0
c        ..........find largest column sum..........
      if (k .le. nopivk) go to 50
         do 30 j=k,n
              if (sigma .ge. sum(j))  go to 30
              sigma = sum(j)
              jbar = j
   30    continue
c
         if (jbar .eq. 0)  go to 220
         if (jbar .eq. k)  go to 50
c        ..........column interchange..........
         i = ipivot(k)
         ipivot(k) = ipivot(jbar)
         ipivot(jbar) = i
         sum(jbar) = sum(k)
         sum(k) = sigma
c
         do 40 i=1,m
              sigma = qr(i,k)
              qr(i,k) = qr(i,jbar)
              qr(i,jbar) = sigma
   40    continue
c        ..........end of column interchange..........
   50    continue
c        ..........  second inner product  ..........
         qrkmax = zero
c
         do 60 i=k,m
              if (dabs( qr(i,k) ) .gt. qrkmax)  qrkmax = dabs( qr(i,k) )
   60    continue
c
         if (qrkmax .lt. ufeta)  go to 210
         alphak = v2norm(mk1, qr(k,k)) / qrkmax
         sigma = alphak**2
c
c        ..........  end second inner product  ..........
         qrkk = qr(k,k)
         if (qrkk .ge. zero)  alphak = -alphak
         alpha(k) = alphak * qrkmax
         beta = qrkmax * dsqrt(sigma - (qrkk*alphak/qrkmax) )
         qr(k,k) = qrkk - alpha(k)
         do 65 i=k,m
   65         qr(i,k) =  qr(i,k) / beta
         k1 = k + 1
         if (k1 .gt. n) go to 120
c
         do 110 j = k1, n
              temp = -dotprd(mk1, qr(k,k), qr(k,j))
c
c             ***  set qr(i,j) = qr(i,j) + temp*qr(i,k), i = k,...,m.
c
              call vaxpy(mk1, qr(k,j), temp, qr(k,k), qr(k,j))
c
              if (k1 .gt. m) go to 110
              sumj = sum(j)
              if (sumj .lt. ufeta) go to 110
              temp = dabs(qr(k,j)/sumj)
              if (temp .lt. rktol1) go to 110
              if (temp .ge. p99) go to 90
                   sum(j) = sumj * dsqrt(one - temp**2)
                   go to 110
   90         sum(j) = v2norm(m-k, qr(k1,j))
  110    continue
c        ..........end of k-th householder transformation..........
  120 continue
c
      go to 999
c     ..........error exit on k-th transformation..........
  210 ierr = -k
      go to 230
c     ..........no non-zero acceptable pivot found..........
  220 ierr = k
  230 do 240 i = k, n
         alpha(i) = zero
         if (i .gt. k) call vscopy(i-k, qr(k,i), zero)
 240     continue
c     ..........return to caller..........
  999 return
c     ..........last card of qrfact..........
      end
      double precision function reldst(p, d, x, x0)
c
c  ***  compute and return relative difference between x and x0  ***
c  ***  nl2sol version 2.2  ***
c
      integer p
      double precision d(p), x(p), x0(p)
c/+
      double precision dabs
c/
      integer i
      double precision emax, t, xmax, zero
c/6
      data zero/0.d+0/
c/7
c     parameter (zero=0.d+0)
c/
c
      emax = zero
      xmax = zero
      do 10 i = 1, p
         t = dabs(d(i) * (x(i) - x0(i)))
         if (emax .lt. t) emax = t
         t = d(i) * (dabs(x(i)) + dabs(x0(i)))
         if (xmax .lt. t) xmax = t
 10      continue
      reldst = zero
      if (xmax .gt. zero) reldst = emax / xmax
 999  return
c  ***  last card of reldst follows  ***
      end
      subroutine rptmul(func, ipivot, j, nn, p, rd, x, y, z)
c
c  ***  func = 1... set  y = rmat * (perm**t) * x.
c  ***  func = 2... set  y = perm * (rmat**t) * rmat * (perm**t) * x.
c  ***  func = 3... set  y = perm * (rmat**t) x.
c
c
c  ***  perm = matrix whose i-th col. is the ipivot(i)-th unit vector.
c  ***  rmat is the upper triangular matrix whose strict upper triangle
c  ***       is stored in  j  and whose diagonal is stored in rd.
c  ***  z is a scratch vector.
c  ***  x and y may share storage.
c
      integer func, nn, p
      integer ipivot(p)
      double precision j(nn,p), rd(p), x(p), y(p), z(p)
c
c  ***  local variables  ***
c
      integer i, im1, k, km1
      double precision zk
c
c  ***  external function  ***
c
      external dotprd
      double precision dotprd
c
c-----------------------------------------------------------------------
c
      if (func .gt. 2) go to 50
c
c  ***  first set  z = (perm**t) * x  ***
c
      do 10 i = 1, p
         k = ipivot(i)
         z(i) = x(k)
 10      continue
c
c  ***  now set  y = rmat * z  ***
c
      y(1) = z(1) * rd(1)
      if (p .le. 1) go to 40
      do 30 k = 2, p
         km1 = k - 1
         zk = z(k)
         do 20 i = 1, km1
 20           y(i) = y(i) + j(i,k)*zk
         y(k) = zk*rd(k)
 30      continue
c
 40   if (func .le. 1) go to 999
      go to 70
c
 50   do 60 i = 1, p
 60      y(i) = x(i)
c
c  ***  set  z = (rmat**t) * y  ***
c
 70   z(1) = y(1) * rd(1)
      if (p .eq. 1) go to 90
      do 80 i = 2, p
         im1 = i - 1
         z(i) = y(i)*rd(i) + dotprd(im1, j(1,i), y)
 80      continue
c
c  ***  now set  y = perm * z  ***
c
 90   do 100 i = 1, p
         k = ipivot(i)
         y(k) = z(i)
 100     continue
c
 999  return
c  ***  last card of rptmul follows  ***
      end
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
      data half/0.5d+0/, one/1.d+0/, zero/0.d+0/
c/7
c     parameter (half=0.5d+0, one=1.d+0, zero=0.d+0)
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
      subroutine slvmul(p, y, s, x)
c
c  ***  set  y = s * x,  s = p x p symmetric matrix.  ***
c  ***  lower triangle of  s  stored rowwise.         ***
c
c  ***  parameter declarations  ***
c
      integer p
      double precision s(1), x(p), y(p)
c     dimension s(p*(p+1)/2)
c
c  ***  local variables  ***
c
      integer i, im1, j, k
      double precision xi
c
c  ***  no intrinsic functions  ***
c
c  ***  external function  ***
c
      external dotprd
      double precision dotprd
c
c-----------------------------------------------------------------------
c
      j = 1
      do 10 i = 1, p
         y(i) = dotprd(i, s(j), x)
         j = j + i
 10      continue
c
      if (p .le. 1) go to 999
      j = 1
      do 40 i = 2, p
         xi = x(i)
         im1 = i - 1
         j = j + 1
         do 30 k = 1, im1
              y(k) = y(k) + s(j)*xi
              j = j + 1
 30           continue
 40      continue
c
 999  return
c  ***  last card of slvmul follows  ***
      end
      logical function stopx(idummy)
c     *****parameters...
      integer idummy
c
c     ..................................................................
c
c     *****purpose...
c     this function may serve as the stopx (asynchronous interruption)
c     function for the nl2sol (nonlinear least-squares) package at
c     those installations which do not wish to implement a
c     dynamic stopx.
c
c     *****algorithm notes...
c     at installations where the nl2sol system is used
c     interactively, this dummy stopx should be replaced by a
c     function that returns .true. if and only if the interrupt
c     (break) key has been pressed since the last call on stopx.
c
c     ..................................................................
c
      stopx = .false.
      return
      end
      subroutine vaxpy(p, w, a, x, y)
c
c  ***  set w = a*x + y  --  w, x, y = p-vectors, a = scalar  ***
c
      integer p
      double precision a, w(p), x(p), y(p)
c
      integer i
c
      do 10 i = 1, p
 10      w(i) = a*x(i) + y(i)
      return
      end
      subroutine vcopy(p, y, x)
c
c  ***  set y = x, where x and y are p-vectors  ***
c
      integer p
      double precision x(p), y(p)
c
      integer i
c
      do 10 i = 1, p
 10      y(i) = x(i)
      return
      end
      subroutine vscopy(p, y, s)
c
c  ***  set p-vector y to scalar s  ***
c
      integer p
      double precision s, y(p)
c
      integer i
c
      do 10 i = 1, p
 10      y(i) = s
      return
      end
      double precision function v2norm(p, x)
c
c  ***  return the 2-norm of the p-vector x, taking  ***
c  ***  care to avoid the most likely underflows.    ***
c
      integer p
      double precision x(p)
c
      integer i, j
      double precision one, r, scale, sqteta, t, xi, zero
c/+
      double precision dabs, dsqrt
c/
      external rmdcon
      double precision rmdcon
c
c/6
      data one/1.d+0/, zero/0.d+0/
c/7
c     parameter (one=1.d+0, zero=0.d+0)
c     save sqteta
c/
      data sqteta/0.d+0/
c
      if (p .gt. 0) go to 10
         v2norm = zero
         go to 999
 10   do 20 i = 1, p
         if (x(i) .ne. zero) go to 30
 20      continue
      v2norm = zero
      go to 999
c
 30   scale = dabs(x(i))
      if (i .lt. p) go to 40
         v2norm = scale
         go to 999
 40   t = one
      if (sqteta .eq. zero) sqteta = rmdcon(2)
c
c     ***  sqteta is (slightly larger than) the square root of the
c     ***  smallest positive floating point number on the machine.
c     ***  the tests involving sqteta are done to prevent underflows.
c
      j = i + 1
      do 60 i = j, p
         xi = dabs(x(i))
         if (xi .gt. scale) go to 50
              r = xi / scale
              if (r .gt. sqteta) t = t + r*r
              go to 60
 50           r = scale / xi
              if (r .le. sqteta) r = zero
              t = one  +  t * r*r
         scale = xi
 60      continue
c
      v2norm = scale * dsqrt(t)
 999  return
c  ***  last card of v2norm follows  ***
      end
c///////////////////////////////////////////////////////////////////////
c  ***  run nl2sol on various test problems, print summary statistics.
c
c     *****common storage with nltest.
c
      common /testcm/ v, rs, jac, nout, nprob, xscal1, xscal2, is, iv
      common /testch/ name, irc
      integer is(6,50), iv(80), jac, nout, nprob, xscal1, xscal2
      real rs(5,50)
c/6
      real name(2,50)
      integer irc(50)
c/7
c     character name(2,50)*4, irc(50)*1
c/
      double precision v(1736)
c
c
c     ..................................................................
c
c     *****purpose.
c        this main program calls nltest to run nl2sol, the nonlinear
c     least-squares solver of ref. 1, on various test problems.
c
c
c     *****application and usage restrictions.
c     this main driver is intended to check whether the nl2sol
c     (nonlinear least-squares) package was successfully
c     transported to a new machine.
c
c     *****algorithm notes.
c     the test problems used are from references (2), (3), and (4).
c     some additional test problems were suggested by jorge more (pri-
c     vate communication).  calls passing these problems to nltest have
c     been commented out (since there are enough other problems), but
c     not removed, since they may be of interest to other researchers.
c
c     *****functions and subroutines called.
c
c        dfault - establishes the default parameter settings for
c                 iv and v.
c
c        imdcon - imdcon(2) returns i/o unit number on which nltest
c                  writes a summary of each test run.
c
c        ivvset - supplies nondefault values for iv and v.
c
c        nltest - calls nl2sol, the nonlinear least-squares
c                  problem solver.
c
c        today  - supplies date and time (or current version of nl2sol).
c
c     *****references.
c
c     (1). dennis, j.e.. gay, d.m.. and welsch, r.e. (1980),
c          an adaptive nonlinear least-squares algorithm,
c          submitted to acm trans. math. software.
c          under revision.
c
c     (2). gill, p.e.. and murray, w. (1976),algorithms for the
c          solution of the non-linear least-squares problem,
c          npl report nac71,(national physical laboratory,
c          division of numerical analysis and computing,
c          teddington,middlesex,england).
c
c     (3) meyer, r.r. (1970), theoretical and computational aspects
c        of nonlinear regression, pp. 465-486 of nonlinear programming,
c        edited by j.b. rosen, o.l.mangasarian, and k. ritter,
c        academic press, new york.
c
c     (4) brown, k.m. (1969), a quadratically convergent newton-
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
c     *****intrinsic functions.
c/+
      integer mod
      double precision dmax1
c/
c     *****external functions and subroutines.
      external dfault, imdcon, ivvset, nltest, today
      integer imdcon
c
c     *****local variables.
      logical rstart
      integer i, j, k, mxfcsv, mxitsv, pu
c/6
      integer jtyp(2)
      real datime(4)
c/7
c     character datime(4)*4, jtyp(2)*1
c/
c
c/6
      data rstart/.false./, jtyp(1),jtyp(2)/1h ,1h*/
c/7
c     data rstart/.false./, jtyp(1),jtyp(2)/' ','*'/
c/
c
c-----------------------------------------------------------------------
c
c  ***  establish default parameter settings  ***
      call dfault (iv, v)
      nout = imdcon(2)
c
c  ***  non-default parameter settings  ***
c
      call ivvset(iv, v)
      pu = iv(21)
c
      jac = 1
      nprob = 0
      xscal1 = 1
      xscal2 = 3
c
c/6
      call nltest(2,2,1,4hrosn,4hbrok,rstart)
      call nltest(3,3,2,4hheli,4hx   ,rstart)
      call nltest(4,4,3,4hsing,4hular,rstart)
      call nltest(7,4,4,4hwood,4hs   ,rstart)
      xscal2 = 1
      call nltest(3,3,5,4hzang,4hwill,rstart)
      xscal2 = 3
      call nltest(5,3,6,4hengv,4hall ,rstart)
      call nltest(2,2,7,4hbran,4hin  ,rstart)
      xscal2 = 2
      call nltest(3,2,8,4hbeal,4he   ,rstart)
      call nltest(5,4,9,4hcrag,4hg   ,rstart)
      xscal2 = 2
      call nltest(10,3,10,4hbox ,4h    ,rstart)
      mxfcsv = iv(17)
      mxitsv = iv(18)
      iv(17) = 20
      iv(18) = 15
      xscal2 = 1
      call nltest(15,15,11,4hdavi,4hdon1,rstart)
      iv(17) = mxfcsv
      iv(18) = mxitsv
      xscal2 = 3
      call nltest(2,2,12,4hfrds,4htein,rstart)
      xscal2 = 1
      call nltest(31,6,13,4hwats,4hon6 ,rstart)
      call nltest(31,9,14,4hwats,4hon9 ,rstart)
      call nltest(31,12,15,4hwats,4hon12,rstart)
      mxfcsv = iv(17)
      iv(17) = 20
      mxitsv = iv(18)
      iv(18) = 15
      call nltest(31,20,16,4hwats,4hon20,rstart)
      iv(17) = mxfcsv
      iv(18) = mxitsv
      xscal2 = 2
      call nltest(8,8,17,4hcheb,4hqd8 ,rstart)
      xscal2 = 3
      call nltest(20,4,18,4hbrow,4hn   ,rstart)
      call nltest(15,3,19,4hbard,4h    ,rstart)
      xscal2 = 1
      call nltest(10,2,20,4hjenn,4hrich,rstart)
      xscal2 = 3
      call nltest(11,4,21,4hkowa,4hlik ,rstart)
      xscal2 = 1
      call nltest(33,5,22,4hosbo,4hrne1,rstart)
      xscal2 = 2
      call nltest(65,11,23,4hosbo,4hrne2,rstart)
      xscal2 = 3
      call nltest(3,2,24,4hmads,4hen  ,rstart)
      xscal2 = 1
      iv(17) = 400
      iv(18) = 300
      call nltest(16,3,25,4hmeye,4hr   ,rstart)
c/7
c     call nltest(2,2,1,'rosn','brok',rstart)
c     call nltest(3,3,2,'heli','x   ',rstart)
c     call nltest(4,4,3,'sing','ular',rstart)
c     call nltest(7,4,4,'wood','s   ',rstart)
c     xscal2 = 1
c     call nltest(3,3,5,'zang','will',rstart)
c     xscal2 = 3
c     call nltest(5,3,6,'engv','all ',rstart)
c     call nltest(2,2,7,'bran','in  ',rstart)
c     xscal2 = 2
c     call nltest(3,2,8,'beal','e   ',rstart)
c     call nltest(5,4,9,'crag','g   ',rstart)
c     xscal2 = 2
c     call nltest(10,3,10,'box ','    ',rstart)
c     mxfcsv = iv(17)
c     mxitsv = iv(18)
c     iv(17) = 20
c     iv(18) = 15
c     xscal2 = 1
c     call nltest(15,15,11,'davi','don1',rstart)
c     iv(17) = mxfcsv
c     iv(18) = mxitsv
c     xscal2 = 3
c     call nltest(2,2,12,'frds','tein',rstart)
c     xscal2 = 1
c     call nltest(31,6,13,'wats','on6 ',rstart)
c     call nltest(31,9,14,'wats','on9 ',rstart)
c     call nltest(31,12,15,'wats','on12',rstart)
c     mxfcsv = iv(17)
c     iv(17) = 20
c     mxitsv = iv(18)
c     iv(18) = 15
c     call nltest(31,20,16,'wats','on20',rstart)
c     iv(17) = mxfcsv
c     iv(18) = mxitsv
c     xscal2 = 2
c     call nltest(8,8,17,'cheb','qd8 ',rstart)
c     xscal2 = 3
c     call nltest(20,4,18,'brow','n   ',rstart)
c     call nltest(15,3,19,'bard','    ',rstart)
c     xscal2 = 1
c     call nltest(10,2,20,'jenn','rich',rstart)
c     xscal2 = 3
c     call nltest(11,4,21,'kowa','lik ',rstart)
c     xscal2 = 1
c     call nltest(33,5,22,'osbo','rne1',rstart)
c     xscal2 = 2
c     call nltest(65,11,23,'osbo','rne2',rstart)
c     xscal2 = 3
c     call nltest(3,2,24,'mads','en  ',rstart)
c     xscal2 = 1
c     iv(17) = 400
c     iv(18) = 300
c     call nltest(16,3,25,'meye','r   ',rstart)
c/
c  ***  brown5  ***
c     call nltest(5,5,26,4hbrow,4hn5  ,rstart)
c  ***  brown10  ***
c     call nltest(10,10,27,4hbrow,4hn10 ,rstart)
c  ***  brown30  ***
c     call nltest(30,30,28,4hbrow,4hn30 ,rstart)
c  ***  brown40  ***
c     call nltest(40,40,29,4hbrow,4hn40 ,rstart)
c  ***  bard+10 ***
c     call nltest(15,3,30,4hbard,4h+10 ,rstart)
c  ***  kowalik and osborne + 10  ***
c     call nltest(11,4,31,4hkowa,4hl+10,rstart)
c  ***  meyer + 10  ***
c     call nltest(16,3,32,4hmeye,4hr+10,rstart)
c  ***  watson6 + 10  ***
c     call nltest(31,6,33,4hwat6,4h+10 ,rstart)
c  ***  watson9 + 10  ***
c     call nltest(31,9,34,4hwat9,4h+10 ,rstart)
c  ***  watson12 + 10  ***
c     call nltest(31,12,35,4hwat1,4h2+10,rstart)
c  ***  watson20 + 10  ***
c     call nltest(31,20,36,4hwat2,4h0+10,rstart)
c
c  ***  repeat two tests using finite-difference jacobian  ***
c
      jac = 2
      xscal2 = 1
c
      iv(17) = 50
      iv(18) = 40
c/6
      call nltest(2,2,1,4hrosn,4hbrok,rstart)
c/7
c     call nltest(2,2,1,'rosn','brok',rstart)
c/
      v(29) = dmax1(1.0d-7, v(29))
      iv(17) = 30
      iv(18) = 20
c  ***  brown  ***
c/6
      call nltest(20,4,18,4hbrow,4hn   ,rstart)
c/7
c     call nltest(20,4,18,'brow','n   ',rstart)
c/
c
      if (nprob .eq. 0 .or. pu .eq. 0) stop
      call today(datime)
      do 130 k = 1, nprob
         if (mod(k,56) .eq. 1) write(pu, 110) datime, nprob
 110     format(1h1,11x,2a4,2x,2a4,10x,10hsummary of,i4,
     1          22h nl2sol test runs.....,10x,
     2          32h(* = finite-difference jacobian)/
     3          48h0 problem    n   p  niter   nf   ng  iv1  x0scal,5x,
     4          39hfinal f     preldf     nreldf     reldx/)
         j = is(6,k)
         write(pu,120) jtyp(j), name(1,k), name(2,k),
     1                 (is(i,k), i=1,5), irc(k), (rs(i,k), i=1,5)
 120     format(1x,a1,2a4,2i4,i7,2i5,3x,a1,f9.1,e13.3,3e11.3)
 130     continue
c
      stop
c...... last card of nlmain ............................................
      end
      subroutine ivvset(iv, v)
c
c  ***  supply nondefault iv and v values for nlmain  (nl2sol ver. 2.2).
c
      integer iv(24)
      double precision v(100)
c
c     activate the next line to turn off detailed summary printing
c     iv(21) = 0
      return
      end
      subroutine nltest (n, p, nex, title1, title2, rstart)
c
c  ***  call nl2sol, save and print statistics  ***
c
c
      integer n, p, nex
      logical rstart
c/6
      real title1, title2
c/7
c     character*4 title1, title2
c/
c
      common /testcm/ v, rs, jac, nout, nprob, xscal1, xscal2, is, iv
      common /testch/ name, irc
      integer is(6,50), iv(80), jac, nout, nprob, xscal1, xscal2
      real rs(5,50)
c/6
      integer irc(50)
      real name(2,50)
c/7
c     character name(2,50)*4, irc(50)*1
c/
      double precision v(1736)
c
      logical rstrt
      integer i, irun, pu, uip(1)
c/6
      integer alg(2), jtyp(2), rc(10)
      real datime(4)
c/7
c     character*4 datime(4)
c     character*2 alg(2)
c     character*1 jtyp(2), rc(10)
c/
      double precision one, t, urparm(1), x(20), x0scal, zero
c
c     ***  external functions and subroutines  ***
c
      external nl2sno, nl2sol, testr, testj, today, xinit
c
c  ***  iv and v subscripts  ***
c
      integer f, f0, nfcall, nfcov, ngcall, niter, nreduc, preduc,
     1        prunit, reldx
c
c/6
      data f/10/, f0/13/, nfcall/6/, nfcov/40/, ngcall/30/,
     1     ngcov/41/, niter/31/, nreduc/6/, preduc/7/,
     2     prunit/21/, reldx/17/
c/7
c     parameter (f=10, f0=13, nfcall=6, nfcov=40, ngcall=30,
c    1     ngcov=41, niter=31, nreduc=6, preduc=7,
c    2     prunit=21, reldx=17)
c/
c/6
      data one/1.d+0/, zero/0.d+0/
c/7
c     parameter (one=1.d+0, zero=0.d+0)
c/
c/6
      data alg(1),alg(2)/2hol,2hno/, jtyp(1),jtyp(2)/1h ,1h*/
      data rc(1)/1h./, rc(2)/1h+/, rc(3)/1hx/, rc(4)/1hr/, rc(5)/1hb/,
     1     rc(6)/1ha/, rc(7)/1hs/, rc(8)/1hf/, rc(9)/1he/, rc(10)/1hi/
c/7
c     data alg(1),alg(2)/'ol','no'/, jtyp(1),jtyp(2)/' ','*'/
c     data rc(1)/'.'/, rc(2)/'+'/, rc(3)/'x'/, rc(4)/'r'/, rc(5)/'b'/,
c    1     rc(6)/'a'/, rc(7)/'s'/, rc(8)/'f'/, rc(9)/'e'/, rc(10)/'i'/
c/
c
c-----------------------------------------------------------------------
c
      uip(1) = nex
      rstrt = rstart
      if (rstrt) go to 20
         pu = iv(prunit)
         call today(datime)
         if (pu .ne. 0) write(pu,10) alg(jac), title1, title2, datime
 10      format (1h1//11h ***** nl2s,a2,12h on problem ,2a4,6h *****,6x,
     1           2a4,2x,2a4)
c
 20   do 100 irun = xscal1, xscal2
         if (rstrt) go to 40
         iv(1) = 12
         x0scal = 1.0d1 ** (irun-1)
c
c        ***  initialize the solution vector x  ***
         call xinit(p, x, nex)
         do 30 i = 1, p
 30           x(i) = x0scal * x(i)
c
 40      if (jac .eq. 1)
     1             call nl2sol(n,p,x,testr,testj,iv,v,uip,urparm,testr)
         if (jac .eq. 2)
     1             call nl2sno(n,p,x,testr,iv,v,uip,urparm,testr)
         if (.not. rstrt .and. nprob .lt. 50) nprob = nprob + 1
         name(1,nprob) = title1
         name(2,nprob) = title2
         is(1,nprob) = n
         is(2,nprob) = p
         is(3,nprob) = iv(niter)
         is(4,nprob) = iv(nfcall) - iv(nfcov)
         is(5,nprob) = iv(ngcall) - iv(ngcov)
         i = iv(1)
         irc(nprob) = rc(i)
         is(6,nprob) = jac
         rs(1,nprob) = x0scal
         rs(2,nprob) = v(f)
         t = one
         if (v(f0) .gt. zero) t = v(preduc) / v(f0)
         rs(3,nprob) = t
         t = one
         if (v(f0) .gt. zero) t = v(nreduc) / v(f0)
         rs(4,nprob) = t
         rs(5,nprob) = v(reldx)
         rstrt = .false.
         if (nout .eq. 0) go to 100
         if (nprob .eq. 1) write(nout,50) datime
 50      format(1h1,11x,2a4,2x,2a4,10x,24hnl2sol test summary.....,10x,
     1          32h(* = finite-difference jacobian)/
     2          48h0 problem    n   p  niter   nf   ng  iv1  x0scal,5x,
     3          39hfinal f     preldf     nreldf     reldx/)
         write(nout,60) jtyp(jac), title1, title2,
     1                (is(i,nprob),i=1,5),irc(nprob),(rs(i,nprob),i=1,5)
 60      format(1x,a1,2a4,2i4,i7,2i5,3x,a1,f9.1,e13.3,3e11.3)
 100     continue
c
 999  return
c  ***  last card of nltest follows  ***
      end
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
      data twopi/6.283185307179586d+0/, zero/0.d+0/
c/7
c     parameter (twopi=6.283185307179586d+0, zero=0.d+0)
c/
c/6
c/7
c     save expmin, uftolg
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
      data twopi/6.283185307179586d+0/
c/7
c     parameter (twopi=6.283185307179586d+0)
c/
c/6
c/7
c     save expmax, expmin, uftolg
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
      subroutine today(datime)
c
c  ***  supply sumsol version  ***
c
c/6
      real datime(4), dt1, dt2, dt3, dt4
      data dt1,dt2,dt3,dt4/4hnl2s,4hol  ,4hver.,4h2.2 /
c/7
c     character*4 datime(4), dt1, dt2, dt3, dt4
c     data dt1,dt2,dt3,dt4/'nl2s','ol  ','ver.','2.2 '/
c/
c
      datime(1) = dt1
      datime(2) = dt2
      datime(3) = dt3
      datime(4) = dt4
 999  return
c  ***  last line of datime follows  ***
      end
      subroutine xinit(p, x, nex)
c
c     *****parameters...
c
      integer nex, p
      double precision x(p)
c
c     ..................................................................
c
c     *****purpose...
c     this routine initializes the solution vector x according to
c     the initial values for the various test functions given in
c     references (1), (2), and (3).
c     subroutines testr and testj.  (see testr for references.)
c
c     *****parameter description...
c     on input...
c
c        nex is the test problem number.
c
c        p is the number of parameters.
c
c     on output...
c
c        x is the initial guess to the solution.
c
c     *****application and usage restrictions...
c     this routine is called by nltest.
c
c     ..................................................................
c
c     *****local variables...
      integer i
      double precision pp1inv
c     *****intrinsic functions...
c/+
      real float
      double precision dble
c/
      dfloat(ii) = dble(float(ii))
c
      go to (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100,
     1   1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100,
     2   2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 1900, 2100,
     3   2500, 1300, 1400, 1500, 1600),nex
c
c  ***  rosenbrock  ***
 100  x(1) = -1.2d0
      x(2) = 1.0d0
      go to 9999
c  ***  helix  ***
 200  x(1) = -1.0d0
      x(2) = 0.0d0
      x(3) = 0.0d0
      go to 9999
c  *** singular  ***
 300  x(1) = 3.0d0
      x(2) = -1.0d0
      x(3) = 0.0d0
      x(4) = 1.0d0
      go to 9999
c  ***  woods  ***
 400  x(1) = -3.0d0
      x(2) = -1.0d0
      x(3) = -3.0d0
      x(4) = -1.0d0
      go to 9999
c  ***  zangwill  ***
 500  x(1) = 1.0d2
      x(2) = -1.0d0
      x(3) = 2.5d0
      go to 9999
c  ***  engvall  ***
 600  x(1) = 1.0d0
      x(2) = 2.0d0
      x(3) = 0.0d0
      go to 9999
c  *** branin  ***
 700  x(1) = 2.0d0
      x(2) = 0.0d0
      go to 9999
c  ***  beale  ***
 800  x(1) = 1.0d-1
      x(2) = 1.0d-1
      go to 9999
c  *** cragg and levy  ***
 900  x(1) = 1.0d0
      x(2) = 2.0d0
      x(3) = 2.0d0
      x(4) = 2.0d0
      go to 9999
c  ***  box  ***
 1000 x(1) = 0.0d0
      x(2) = 1.0d1
      x(3) = 2.0d1
      go to 9999
c  ***  davidon 1  ***
 1100 do 1101 i = 1,p
 1101    x(i) = 0.0d0
      go to 9999
c  ***  freudenstein and roth  ***
 1200 x(1) = 1.5d1
      x(2) = -2.0d0
      go to 9999
c  ***  watson  ***
 1300 continue
 1400 continue
 1500 continue
 1600 do 1601 i = 1,p
 1601    x(i) = 0.0d0
      go to 9999
c  ***  chebyquad  ***
 1700 pp1inv = 1.0d0/dfloat(p + 1)
      do 1701 i = 1, p
 1701    x(i) = dfloat(i)*pp1inv
      go to 9999
c  *** brown and dennis  ***
 1800 x(1) = 2.5d1
      x(2) = 5.0d0
      x(3) = -5.0d0
      x(4) = -1.0d0
      go to 9999
c  ***  bard  ***
 1900 x(1) = 1.d0
      x(2) = 1.d0
      x(3) = 1.d0
      go to 9999
c  ***  jennrich and sampson  ***
 2000 x(1) = 3.0d-1
      x(2) = 4.0d-1
      go to 9999
c  ***  kowalik and osborne  ***
 2100 x(1) = 2.5d-1
      x(2) = 3.9d-1
      x(3) = 4.15d-1
      x(4) = 3.9d-1
      go to 9999
c  ***  osborne 1  ***
 2200 x(1) = 5.0d-1
      x(2) = 1.5d0
      x(3) = -1.0d0
      x(4) = 1.0d-2
      x(5) = 2.0d-2
      go to 9999
c  ***  osborne 2  ***
 2300 x(1) = 1.3d0
      x(2) = 6.5d-1
      x(3) = 6.5d-1
      x(4) = 7.0d-1
      x(5) = 6.0d-1
      x(6) = 3.0d0
      x(7) = 5.0d0
      x(8) = 7.0d0
      x(9) = 2.0d0
      x(10) = 4.5d0
      x(11) = 5.5d0
      go to 9999
c  ***  madsen  ***
 2400 x(1) = 3.0d0
      x(2) = 1.0d0
      go to 9999
c  ***  meyer  **
 2500 x(1) = 2.0d-2
      x(2) = 4.0d3
      x(3) = 2.5d2
      go to 9999
c  ***  brown  ***
 2600 continue
 2700 continue
 2800 continue
 2900 do 2901 i = 1, p
 2901    x(i) = 5.d-1
      go to 9999
c
c
 9999 return
      end
