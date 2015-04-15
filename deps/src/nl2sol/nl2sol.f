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
c      data nfcall/6/, nfgcal/7/, toobig/2/
c/7
      parameter (nfcall=6, nfgcal=7, toobig=2)
c/
c
c  ***  v subscript values  ***
c
c/6
c      data d/27/, j/33/, r/50/
c/7
      parameter (d=27, j=33, r=50)
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
