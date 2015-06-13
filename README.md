# NL2SOL: Non-linear least squares optimization

This is the original netlib version of NL2SOL, a non-linear,
least-squares optimization program.  It is detailed in two
Transactions on Mathematical Software (TOMS) papers.  They are:

J.E. Dennis, D.M. Gay, R.E. Welsch, "An Adaptive Nonlinear
Least-Squares Algorithm", ACM Transactions on Mathematical Software
(TOMS), Volume 7 Issue 3, Sept. 1981, pp 348-368, ACM New York, NY, USA
[see here](http://dl.acm.org/citation.cfm?id=355965&CFID=660003329&CFTOKEN=25049918)

J.E. Dennis, D.M. Gay, R.E. Welsch, "Algorithm 573: NL2SOL—An Adaptive
Nonlinear Least-Squares Algorithm", ACM Transactions on Mathematical
Software (TOMS), Volume 7 Issue 3, Sept. 1981, pp 369-383, ACM New
York, NY, USA [see here](http://dl.acm.org/citation.cfm?id=355966)

Here we use the original NL2SOL Fortran 77 source code which appears
as TOMS ALgorithm 573 (NL2SOL Version 2.2).  The code was downloaded
from netlib and is archived in the deps/src/nl2sol directory as a
single blob in the file named nl2sol.netlib.orig.f

This blob has also been broken up into the individual source files and
commented out the "c/6" code for the "c/7" code, which enables the f77
version.  Also added are cmake files for building the code and running
the tests.  Running the Fortran tests and coverage is manual and is
not part of the installation. (The coverage is a very respectable 87%)
The original fortran test code now lives in a separate subdirectory
(.../deps/src/tests) as well.

Wrapper code has been added using the C interface facilities of Julia.
(ie ccall and cfunction etc), so that nl2sol can be called directly
from Julia.

The runtests.jl in the test directory has many examples of using Julia
to call nl2sol and using Julia functions to calculate the residual and
the jacobian.  These reproduce the tests included in the test suite.
Note that the results from the original Fortran test driver and the
newer Julia based residual and jacobian calculations reproduce the
results that are published in the above papers.  Also included there
is the ability to run the same tests by using
Optim.levenberg_marquardt.

There are many tuning parameters for NL2SOL but they have not been
made visible in the calling signiture.  Likewise, NL2SNO, which
will use finite differences to calculate the jacobian has also not
been made visible.  If there is a demand for it, it can certainly
be added.

As an optimization solution, this would compete most directly with the
levenberg\_marquardt from the Optim module.  It differs from that
algorithm in that NL2SOL is a quasi-Newton method (_not_ BFGS but
rather DFP for those who care).  Because of that you would expect
NL2SOL to perform better on those models that have large(r) residuals
at the optimum.  It will also generally perform better if the starting
guess is far from the optimim point.

## Limitations

  * Only supported in Julia 0.4 (and greater, someday)

  * Only a linux version, compiled on Ubuntu 14.04 is currently available.

  * NL2sol does nasty pointer tricks that are liable to make programs
brittle, so don't do a "using NL2sol" in code you care about.

  * nl2sno, which calculates the jacobian by finite differences, has not
been exported.

  * nl2itr, which uses "reverse communication" to request residual and jacobian
updates, has not been exported.

  * Many, many tuning parameters have not been exported.

  * NL2sol uses a different convergence testing strategy than Optim.levenberg_marquardt.
This makes doing apples to apples comparisons challenging.


## Example Usage

Here is a simple and complete example of using NL2SOl.


    using NL2sol

    function rosenbrock_res(x, r)
        r[1] = 10.0 * (x[2] - x[1]^2 )
        r[2] = 1.0 - x[1]
        return r
    end

    function rosenbrock_jac(x, jac)
        jac[1, 1] = -20.0 * x[1]
        jac[1, 2] =  10.0
        jac[2, 1] =  -1.0
        jac[2, 2] =   0.0
       return jac
    end

    function main()
        println("NL2SOL on Rosenbrock")
        result = nl2sol(rosenbrock_res, rosenbrock_jac, [-1.2, 1.0], 2)
        println(result)
    end

    main()


Note that we let the Julia wrapper for nl2sol to allocate the memory
for both the residual and the jacobian.

nl2sol can print detailed iteration summaries.  This is turned on by
setting the keyword parameter quiet to false, ie

        result = nl2sol(rosenbrock_res, rosenbrock_jac, [-1.2, 1.0], 2; quiet=false)
