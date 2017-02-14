# NL2sol.jl: Non-linear least squares optimization

NL2sol.jl solves the non-linear least squares problem.  That is, it
finds an x that minimizes $ \sum_{i=1}^{n} {{r}_{i}}^{2}(x) $ where x
is a vector of size p.  It returns a struct of type
Optim.MultivariateOptimizationResults that contains the relevant info
(see the Julia Optim module docs for further info).  It does this by
wrapping the FORTRAN version of the code.

The residual and the jacobian functions are expected to take args that
have been preallocated for those values.  These arrays are actually
allocated in the Julia function nl2sol before passing to the Fortran
subroutine nl2sol.

NOTE: NL2sol.jl does some pointer tricks in order to interface with the
      FORTRAN code, so tread lightly if you'd like to modify the code.

## EXAMPLE USAGE NL2sol.nl2sol

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
        result = nl2sol(rosenbrock_res, rosenbrock_jac, [-1.2, 1.0], 2; quiet=true)
        println(result)
    end

    main()

Alternatively, if you do not have or do not want to write a jacobian, you can
use nl2sno, which uses a finite difference approximation to the jacobian.  Even
if you do have the jacobian available, nl2sno can be used to check for its 
correctness.  In this case, you must provide the iv and v arrays (see below).
A complete example would look like:

## EXAMPLE USAGE NL2sol.nl2sno

    using NL2sol

    function rosenbrock_res(x, r)
        r[1] = 10.0 * (x[2] - x[1]^2 )
        r[2] = 1.0 - x[1]
        return r
    end

    function main()
        println("NL2nso on Rosenbrock")
        iv, v = nl2_set_defaults(2, 2)
        result = nl2sno(rosenbrock_res, [-1.2, 1.0], 2, iv, v)
        println(result)
    end

    main()


## Background

The wrapped Fortran code is the original netlib version of NL2SOL, a non-linear,
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
the tests.  Running the Fortran tests and coverage is manual and
not part of the installation. (The coverage is a very respectable 87%)
The original fortran test code now lives in a separate subdirectory
(.../deps/src/tests) as well.  To learn how to build and run the tests
with coverage, see NL2sol.jl/deps/src/CMakeLists.txt

Wrapper code has been added using the C interface facilities of Julia.
(ie ccall and cfunction etc), so that nl2sol can be called directly
from Julia.

The runtests.jl in the test directory has many examples of using Julia
to call nl2sol and using Julia functions to calculate the residual and
the jacobian.  Also included there is the ability to run the same tests by
using Optim.levenberg_marquardt.

There are two calling signitures for nl2sol.  One is the simplified
version used above an its complete version is given by:

    function nl2sol(res::Function, jac::Function, init_x, n; 
                    maxIter=df_maxIter, maxFuncCall=df_maxFuncCall, 
                    tolX=df_tolX,  tolAbsFunc=df_tolAbsFunc,
                    tolRelFunc=df_tolRelFunc, quiet=true)

The required arguments are the function that calculates the residual
vector, the funtion that calculates the jacobian, the initial starting
guess for the non-linear parameters, and the number of 'measurements'
that we are fitting (this is also the length of the residual vector
returned by the res::Function).  The optional arguments control
some (but not all) of the convergence criteria.

The alternative version requires that you first call a function to set
the defaults and that function returns an integer and a real array
which must then be passed to nl2sol.  A calling sequence would look
like

    iv, v = set_defaults(n, p)
    ## change default values inside of iv, v
    results = nl2sol(res, jac, init_x, n, iv, v)

The advantage of this form is that all of the control and tuning
parameters of NL2sol are available by changing some of the values in
the iv and/or v arrays.  Also available are more status values in
these arrays. They are well documented in the 'program paper' above.

As an optimization solution, this would compete most directly with the
levenberg\_marquardt from the Optim module.  It differs from that
algorithm in that NL2SOL is a quasi-Newton method (_not_ BFGS but
rather DFP for those who care).  Because of that you would expect
NL2SOL to perform better on those models that have large(r) residuals
at the optimum.  It will also generally perform better if the starting
guess is far from the optimim point.

## Limitations

  * Only supported in Julia 0.6 (dev).

  * Only a linux version, compiled on Ubuntu 16.04 is currently available. But the cmake
build scripts should work on other Linux machines. Windows might be a challenge.

  * nl2itr, which uses "reverse communication" to request residual and jacobian
updates, has not been exported.

  * NL2sol uses a different convergence testing strategy than Optim.levenberg_marquardt.
This makes doing apples to apples comparisons challenging.

Note that we let the Julia wrapper for nl2sol allocate the memory for
both the residual and the jacobian.

nl2sol can print detailed iteration summaries.  This is turned on by
setting the keyword parameter quiet to false, ie

        result = nl2sol(rosenbrock_res, rosenbrock_jac, [-1.2, 1.0], 2; quiet=false)

