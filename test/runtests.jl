# levenberg_marquardt is actually pretty simple _but_
# LsqFit pulls in a bunch of stuff which is broken
# right now
use_levenberg = "use_levenberg" in ARGS
use_levenberg && using LsqFit.levenberg_marquardt

# the finite difference derivatives do nasty things...
use_nl2sno = "use_nl2sno" in ARGS

# use the installed NL2sol instead of local dev one
use_installed = "use_installed" in ARGS
if use_installed
    using NL2sol
else
    ENV["NL2SOL_LIBPATH"] = "../deps/usr/lib"
    include("../src/NL2sol.jl")
    import Main.NL2sol: nl2sol, nl2sno, nl2_set_defaults, nl2_reset_defaults!,
                        PRUNIT, MXITER, MXFCAL, FUNCT0, NREDUC, NFCALL, NFCOV,
                        NGCALL, NGCOV, FUNCT, RELDX
end

using Base.Test
using Formatting

# Often DataFrames is borked in the dev stream
use_DataFrames = "use_DataFrames" in ARGS
use_DataFrames && using DataFrames

# convert the optim multivariate results to a local version that
# can be pretty printed.
function convert_results(r)
    nl_results = NL2sol.NL2OptimizationResults(string(r.method), r.initial_x,
                                               r.minimizer, r.minimum, r.iterations,
                                               r.iteration_converged, r.x_converged, r.x_tol,
                                               r.f_converged, r.f_tol, r.g_converged, r.g_tol,
                                               r.f_calls, r.g_calls)
    return nl_results
end
                                        
# These are the problems we will run and the starting guesses for
# the optimal solution.  We allocate r and j here only for levenberg_marquardt
# (We need them for the closures).  The Julia function nl2sol will allocate 
# these before calling the fortran version
const problems = Dict( 
     # problem      allocate r  allocate j     init X   rescale init X
    "rosenbrock" => (zeros(2), zeros(2, 2), [-1.2, 1.0], 3),
    "helix" => (zeros(3), zeros(3, 3), [-1.0e0, 0.0e0, 0.0e0], 3),
    "singular" => (zeros(4), zeros(4, 4), [3.0, -1.0, 0.0, 1.0], 3),
    "woods" => (zeros(7), zeros(7, 4), [-3.0, -1.0, -3.0, -1.0], 3),
    "zangwill" => (zeros(3), zeros(3, 3), [100., -1., 2.5], 1),
    "engvall" => (zeros(5), zeros(5, 3), [1., 2., 0.], 3),
    "branin" => (zeros(2), zeros(2, 2), [2., 0.], 3),
    "beale" => (zeros(3), zeros(3, 2), [0.1, 0.1], 2),
    "cragg_levy" => (zeros(5), zeros(5, 4), [1.,2.,2.,2.], 2),
    "box" => (zeros(10), zeros(10, 3), [0.0, 1.0, 2.0], 2),
    "davidon" => (zeros(15), zeros(15, 15), zeros(15), 1),
    "freudenstein_roth" => (zeros(2), zeros(2, 2), [15.0, -2.0], 3),
    "watson6" => (zeros(31), zeros(31, 6), zeros(6), 1),
    "watson9" => (zeros(31), zeros(31, 9), zeros(9), 1),
    "watson12" => (zeros(31), zeros(31, 12), zeros(12), 1),
    "watson20" => (zeros(31), zeros(31, 20), zeros(20), 1),
    "chebyquad" => (zeros(8), zeros(8, 8), collect(1.:8)/9., 2),
    "brown_dennis" => (zeros(20), zeros(20, 4), [25., 5., -5., -1.], 3),
    "bard" => (zeros(15), zeros(15, 3), [1., 1., 1.], 3),
    "jennrich_sampson" => (zeros(10), zeros(10, 2), [0.3, 0.4], 1),
    "kowalik_osborne" => (zeros(11), zeros(11, 4), [.25, .39, .415, .39], 3),
    "osborne1" => (zeros(33), zeros(33, 5), [.5, 1.5, -1.0, .01, .02], 1),
    "osborne2" => (zeros(65), zeros(65, 11), [1.3, .65, .65, .7, .6, 3., 5., 
                                              7., 2., 4.5, 5.5], 2),
    "madsen" => (zeros(3), zeros(3, 2), [3., 1.], 3),
    "meyer" => (zeros(16), zeros(16, 3), [0.02, 4.e3, 2.5e2], 1),
    "brown5" => (zeros(5), zeros(5, 5), 0.5*ones(5), 1),
    "brown10" => (zeros(10), zeros(10, 10), 0.5*ones(10), 1),
    "brown30" => (zeros(30), zeros(30, 30), 0.5*ones(30), 1),
    "brown40" => (zeros(40), zeros(40, 40), 0.5*ones(40), 1)
)


lmrj = Dict()

# We need to wrap these functions for levenberg_marquardt
for (k, v) in problems
    res = Symbol(string(k, "_res"))
    lmres = Symbol(string("wr_", k, "_res"))
    jac = Symbol(string(k, "_jac"))
    lmjac = Symbol(string("wr_", k, "_jac"))
    @eval begin
        function ($lmres)(x)
            ($res)(x, $(v[1]))
        end
        function ($lmjac)(x)
            ($jac)(x, $(v[2]))
        end
    end
    lmrj[k] = (eval(lmres), eval(lmjac))
end


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

rrn = 1e-5 * rand(2)
rjn = 1e-5 * rand(2,2)

function helix_res(x, r)
    theta = atan2(x[2], x[1]) / 2pi

    if x[1] <= 0.0 && x[2] <= 0.0
        theta = theta + 1.0
    end

    r[1] = 10.0 * (x[3] - 10.0 * theta )
    r[2] = 10.0 * (sqrt(x[1]^2 + x[2]^2) - 1.0)
    r[3] = x[3]
    return r
end

function helix_jac(x, jac)
    t = x[1]^2 + x[2]^2
    ti = 100.0 / (2pi * t )
    jac[1, 1] = ti * x[2]
    t = 10.0 / sqrt(t)
    jac[2, 1] = x[1] * t
    jac[3, 1] = 0.0
    jac[1, 2] = -ti * x[1]
    jac[2, 2] = x[2] * t
    jac[3, 2] = 0.0
    jac[1, 3] = 10.0
    jac[2, 3] = 0.0
    jac[3, 3] = 1.0
    return jac
end

function singular_res(x, r)
    r[1] = x[1] + 10.0 * x[2]
    r[2] = sqrt(5.0) * (x[3] - x[4])
    r[3] = (x[2] - 2x[3])^2
    r[4] = sqrt(10.0) * (x[1] - x[4])^2
    return r
end

function singular_jac(x, jac)
    #jac[:] = 0.0
    for k in 1:4
        for i in 1:4
            jac[i, k] = 0.0
        end
    end
    jac[1, 1] = 1.0
    jac[1, 2] = 10.0
    jac[2, 3] = sqrt(5.0)
    jac[2, 4] = -jac[2, 3]
    jac[3, 2] = 2.0 * (x[2] - 2.0 * x[3])
    jac[3, 3] = -2.0 * jac[3, 2]
    jac[4, 1] = sqrt(40.0) * (x[1] - x[4])
    jac[4, 4] = -jac[4, 1]
    return jac
end

function woods_res(x, r)
    r[1] = 10.0 * (x[2] - x[1]^2 )
    r[2] = 1.0 - x[1]
    r[3] = sqrt(90.0) * (x[4] - x[3]^2 )
    r[4] = 1.0 - x[3]
    r[5] = sqrt(9.9) * (x[2] + x[4] - 2.0)
    t = sqrt(0.2)
    r[6] = t * (x[2] - 1.0)
    r[7] = t * (x[4] - 1.0)
    return r
end

function woods_jac(x, jac)
    #jac[:] = 0.0
    for k in 1:4
        for i in 1:7
            jac[i, k] = 0.0
        end
    end
    jac[1, 1] = -20x[1]
    jac[1, 2] = 10.0
    jac[2, 1] = -1.0
    jac[3, 4] = sqrt(90.0)
    jac[3, 3] = -2x[3] * jac[3, 4]
    jac[4, 3] = -1.0
    jac[5, 2] = sqrt(9.9)
    jac[5, 4] = jac[5, 2]
    jac[6, 2] = sqrt(0.2)
    jac[7, 4] = jac[6, 2]
    return jac
end

function davidon_res(x, r)
    for i = 1:14
        r1 = 0.0
        ti = Float64(i)
        t = 1.0
        for j = 1:15
            r1 += t * x[j]
            t *= ti
        end
        r[i] = r1
    end
    r[15] = x[1] - 1.0
    return r
end

function davidon_jac(x, jac)
    for i = 1:14
        ti = Float64(i)
        t = 1.0
        for k = 1:15
            jac[i, k] = t
            t *= ti
        end
    end
    jac[15, 1] = 1.0
    for k = 2:15
        jac[15, k] = 0.0
    end
    return jac
end


function zangwill_res(x, r)
    r[1] = x[1] - x[2] + x[3]
    r[2] = -x[1] + x[2] + x[3]
    r[3] = x[1] + x[2] - x[3]
    return r
end


function zangwill_jac(x, jac)
    #jac[:] = 0.0
    for k = 1:3
        for i = 1:3
            jac[i, k] = 1.0
        end
    end
    jac[1, 2] = -1.0
    jac[2, 1] = -1.0
    jac[3, 3] = -1.0
    return jac
end


function engvall_res(x, r)
    r[1] = x[1]^2 + x[2]^2 + x[3]^2 - 1.0
    r[2] = x[1]^2 + x[2]^2 + (x[3] - 2.0)^2 - 1.0
    r[3] = x[1] + x[2] + x[3] - 1.0
    r[4] = x[1] + x[2] - x[3] + 1.0
    r[5] = x[1]^3 + 3.0 * x[2] ^ 2 + (5.0 * x[3] - x[1] + 1.0)^2 - 36.0
    return r
end


function engvall_jac(x, jac)
    jac[1, 1] = 2.0 * x[1]
    jac[1, 2] = 2.0 * x[2]
    jac[1, 3] = 2.0 * x[3]
    jac[2, 1] = jac[1, 1]
    jac[2, 2] = jac[1, 2]
    jac[2, 3] = 2.0 * (x[3] - 2.0)
    jac[3, 1] = 1.0
    jac[3, 2] = 1.0
    jac[3, 3] = 1.0
    jac[4, 1] = 1.0
    jac[4, 2] = 1.0
    jac[4, 3] = -1.0
    t = 2.0 * (5.0 * x[3] - x[1] + 1.0)
    jac[5, 1] = 3.0 * x[1] ^ 2 - t
    jac[5, 2] = 6.0 * x[2]
    jac[5, 3] = 5.0 * t
    return jac
end


function branin_res(x, r)
    r[1] = 4.0 * (x[1] + x[2])
    r[2] = r[1] + (x[1] - x[2]) * ((x[1] - 2.0) ^ 2 + x[2] ^ 2 - 1.0)
    return r
end


function branin_jac(x, jac)
    jac[1, 1] = 4.0
    jac[1, 2] = 4.0
    jac[2, 1] = 3.0 + (x[1] - 2.0) * (3.0 * x[1] - 2.0 * x[2] - 2.0) + x[2] * x[2]
    jac[2, 2] = 1.0 + 2.0 * (2.0*x[1] - x[2] * x[2]) - (x[1] - x[2])^2
    return jac
end

function beale_res(x, r)
    r[1] = 1.5 - x[1] * (1.0 - x[2])
    r[2] = 2.25 - x[1] * (1.0 - x[2]^2)
    r[3] = 2.625 - x[1] * (1.0 -  x[2]^3)
    return r
end

function beale_jac(x, jac)
    jac[1, 1] = x[2] - 1.0
    jac[1, 2] = x[1]
    jac[2, 1] = x[2] ^ 2 - 1.0
    jac[2, 2] = 2.0 * x[1] * x[2]
    jac[3, 1] = x[2] ^ 3 - 1.0
    jac[3, 2] = 3.0 * x[1] * (x[2] ^ 2)
    return jac
end

function cragg_levy_res(x, r)
    r[1] = (exp(x[1]) - x[2])^2
    r[2] = 10.0 * (x[2] - x[3])^3
    r[3] = (sin(x[3] - x[4]) / cos(x[3] - x[4]) )^2
    r[4] = x[1]^4
    r[5] = x[4] - 1.0
    return r
end

function cragg_levy_jac(x, jac)
    for i = 1:5
        for k = 1:4
            jac[i, k] = 0.0
        end
    end
    t = exp(x[1])
    jac[1, 2] = -2.0 * (t - x[2])
    jac[1, 1] = -t * jac[1, 2]
    jac[2, 2] = 30.0 * (x[2] - x[3])^2
    jac[2, 3] = -jac[2, 2]
    jac[3, 3] = 2.0*sin(x[3] - x[4]) / (cos(x[3] - x[4]))^3
    jac[3, 4] = -jac[3, 3]
    jac[4, 1] = 4.0*x[1]^3
    jac[5, 4] = 1.0
    return jac
end

# It turns out that this works OK for starting points 10*x_init and
# 100*x_init.  Interesting, the Fortran code does _not_ run the case
# 1*x_init.   It looks like it converges to a different minimum from
# that starting point.
const expmax = 1.999 * log(sqrt(0.999*realmax()))
const expmin = 1.999 * log(sqrt(1.001*realmin()))

function box_res(x, r)
    if -expmax > minimum(x)
        throw("Error in box residual")
    end
    for i = 1:10
        ti = -0.1 * Float64(i)
        t1 = ti * x[1]
        t1 > expmin ? e1 = exp(t1) : e1 = 0.0
        t2 = ti*x[2]
        t2 > expmin ? e2 = exp(t2) : e2 = 0.0
        r[i] = (e1 - e2) - x[3]*(exp(ti) - exp(10.0 * ti))
    end
    return r
end

function box_jac(x, jac)
    for i = 1:10
        ti = -0.1 * Float64(i)
        t = x[1] * ti
        t > expmin ? e = exp(t) : e = 0.0
        jac[i, 1] = ti * e
        t = x[2] * ti
        t > expmin ? e = exp(t) : e = 0.0
        jac[i, 2] = -ti * e
        jac[i, 3] = exp(10.0 * ti) - exp(ti)
    end
    return jac
end

function freudenstein_roth_res(x, r)
    r[1] = -1.3e1 + x[1] - 2.0e0*x[2] + 5.0e0*x[2]^2 - x[2]^3
    r[2] = -2.9e1 + x[1] - 1.4e1*x[2] + x[2]^2 + x[2]^3
    return r
end



function freudenstein_roth_jac(x, j)
    j[1,1] = 1.e0
    j[1,2] = -2.e0 + x[2]*(1.e1 - 3.e0*x[2])
    j[2,1] = 1.e0
    j[2,2] = -1.4e1 + x[2]*(2.e0 + 3.e0*x[2])
    return j
end

function watson_res(x, r, p)
    # for watson nex = 13|14|15|16
    for i = 1:29
        ti = Float64(i) / 2.9e1
        r1 = 0.0e0
        r2 = x[1]
        t = 1.0e0
        for j = 2:p
              r1 = r1 + Float64(j-1) * t * x[j]
              t = t * ti
              r2 = r2 + t * x[j]
        end
        r[i] = r1 - r2 * r2 - 1.0e0
    end
    r[30] = x[1]
    r[31] = x[2] - x[1] ^ 2 - 1.0e0
    return r
end

watson6_res(x, r) = watson_res(x, r, 6)
watson9_res(x, r) = watson_res(x, r, 9)
watson12_res(x, r) = watson_res(x, r, 12)
watson20_res(x, r) = watson_res(x, r, 20)

function watson_jac(x, j, p)
    for i = 1:29
        ti = Float64(i) / 2.9e1
        r2 = x[1]
        t = 1.e0
        for k = 2:p
             t = t*ti
             r2 = r2 + t*x[k]
        end
        r2 = -2.e0 * r2
        j[i, 1] = r2
        t = 1.e0
        r2 = ti * r2
        for k = 2:p
            j[i, k] = t * (Float64(k-1) + r2)
            t = t * ti
        end
    end

    for i = 30:31
        for k = 2:p
            j[i, k] = 0.e0
        end
    end
    j[30, 1] = 1.e0
    j[31, 1] = -2.e0*x[1]
    j[31, 2] = 1.e0
    return j
end

watson6_jac(x, j) = watson_jac(x, j, 6)
watson9_jac(x, j) = watson_jac(x, j, 9)
watson12_jac(x, j) = watson_jac(x, j, 12)
# the Fortran driver and test code places a 15 iteration limit on 
# watson20 so we will do so as well to match.
watson20_jac(x, j) = watson_jac(x, j, 20)

function chebyquad_res(x, r)
    n = length(r)
    for i = 1:n
       r[i] = 0.0e0
    end
    for j = 1:n
        tim1 = 1.0e0
        ti = 2.0e0 * x[j] - 1.0e0
        z = ti + ti
        for i = 1:n
            r[i] = r[i] + ti
            tip1 = z * ti - tim1
            tim1 = ti
            ti = tip1
        end
    end
    floatn = Float64(n)
    for i = 1:n
        ti = 0.0e0
        if mod(i,2) == 0
            ti = -1.0e0 / Float64(i*i - 1)
        end
        r[i] = ti - r[i]/floatn
    end
    return r
end

function brown_dennis_res(x, r)
    n = length(r)
    for i = 1:n
        ti = 0.2e0 * Float64(i)
        r[i] = (x[1] + x[2] * ti - exp(ti)) ^ 2 +
               (x[3] + x[4] * sin(ti) - cos(ti))^2
    end
    return r
end

const ybard = [1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1,
               3.5e-1, 3.9e-1, 3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1,
               1.34e0, 2.10e0, 4.39e0]

function bard_res(x, r)
    for i = 1:15
        u = Float64(i)
        v = 1.6e1 - u
        w = min(u, v)
        r[i] = ybard[i] - (x[1] + u / (x[2] * v + x[3] * w))
    end
    return r
end

function bard_jac(x, j)
    for i = 1:15
        j[i, 1] = -1.e0
        u = Float64(i)
        v = 1.6e1 - u
        w = min(u, v)
        t = u/(x[2]*v + x[3]*w)^2
        j[i, 2] = v*t
        j[i, 3] = w*t
    end
    return j
end

function jennrich_sampson_res(x, r)
    for i = 1:10
        ti = Float64(i)
        r[i] = 2.0e0 + 2.0e0*ti - (exp(ti*x[1]) + exp(ti*x[2]))
    end
    return r
end

function jennrich_sampson_jac(x, j)
    for i = 1:10
         ti = Float64(i)
         j[i,1] = -ti * exp(ti*x[1])
         j[i,2] = -ti * exp(ti*x[2])
    end
    return j
end

const ykow = [1.957e-1, 1.947e-1, 1.735e-1, 1.600e-1, 8.44e-2, 6.27e-2,
              4.56e-2,  3.42e-2, 3.23e-2, 2.35e-2, 2.46e-2]

const ukow = [4.0e0, 2.0e0, 1.0e0, 5.0e-1, 2.5e-1, 1.67e-1,
              1.25e-1, 1.0e-1, 8.33e-2, 7.14e-2, 6.25e-2]

function kowalik_osborne_res(x, r)
    for i = 1:11
        r[i] = ykow[i] - x[1]*(ukow[i] ^ 2 + x[2] * ukow[i]) / (ukow[i] ^ 2 + x[3]*ukow[i] + x[4])
    end
    return r
end

function kowalik_osborne_jac(x, j)
    for i = 1:11
         t = -1.e0 / (ukow[i]^2 + x[3] * ukow[i] + x[4])
         j[i, 1] = t * (ukow[i]^2 + x[2] * ukow[i])
         j[i, 2] = x[1] * ukow[i]*t
         t = t * j[i, 1] * x[1]
         j[i, 3] = ukow[i] * t
         j[i, 4] = t
    end
    return j
end

const yosb1 = [8.44e-1, 9.08e-1, 9.32e-1, 9.36e-1, 9.25e-1, 9.08e-1,
               8.81e-1, 8.50e-1, 8.18e-1, 7.84e-1, 7.51e-1, 7.18e-1,
               6.85e-1, 6.58e-1, 6.28e-1, 6.03e-1, 5.80e-1, 5.58e-1,
               5.38e-1, 5.22e-1, 5.06e-1, 4.90e-1, 4.78e-1, 4.67e-1,
               4.57e-1, 4.48e-1, 4.38e-1, 4.31e-1, 4.24e-1, 4.20e-1,
               4.14e-1, 4.11e-1, 4.06e-1]

function osborne1_res(x, r)
    for i = 1:33
        ti = 1.0e1 * Float64(1-i)
        r[i] = yosb1[i] - (x[1] + x[2] * exp(x[4] * ti) + x[3] * exp(x[5] * ti))
    end
    return r
end

function osborne1_jac(x, j)
    for i = 1:33
         ti = 1.0e1 * Float64(1-i)
         j[i, 1] = -1.e0
         j[i, 2] = -exp(x[4] * ti)
         j[i, 3] = -exp(x[5] * ti)
         j[i, 4] = ti * x[2] * j[i, 2]
         j[i, 5] = ti * x[3] * j[i, 3]
    end
    return j
end

const yosb2 = [ 1.366e0, 1.191e0, 1.112e0, 1.013e0, 9.91e-1, 8.85e-1,
                8.31e-1, 8.47e-1, 7.86e-1, 7.25e-1, 7.46e-1, 6.79e-1,
                6.08e-1, 6.55e-1, 6.16e-1, 6.06e-1, 6.02e-1, 6.26e-1,
                6.51e-1, 7.24e-1, 6.49e-1, 6.49e-1, 6.94e-1, 6.44e-1,
                6.24e-1, 6.61e-1, 6.12e-1, 5.58e-1, 5.33e-1, 4.95e-1,
                5.00e-1, 4.23e-1, 3.95e-1, 3.75e-1, 3.72e-1, 3.91e-1,
                3.96e-1, 4.05e-1, 4.28e-1, 4.29e-1, 5.23e-1, 5.62e-1,
                6.07e-1, 6.53e-1, 6.72e-1, 7.08e-1, 6.33e-1, 6.68e-1,
                6.45e-1, 6.32e-1, 5.91e-1, 5.59e-1, 5.97e-1, 6.25e-1,
                7.39e-1, 7.10e-1, 7.29e-1, 7.20e-1, 6.36e-1, 5.81e-1,
                4.28e-1, 2.92e-1, 1.62e-1, 9.8e-2, 5.4e-2]

# uftolg is a machine-dependent constant.  it is just slightly
# larger than the log of the smallest positive machine number.
const uftolg = 1.999e0 * log(sqrt(1.001 * 2.2250738585072014e-308))

function osborne2_res(x, r)
    for i = 1:65
         ti = 0.1e0 * Float64(1-i)
         ri = x[1] * exp(x[5] * ti)
         for j = 2:4
              t = 0.e0
              theta = -x[j+4] * (ti + x[j+7]) ^ 2
              if theta > uftolg
                  t = exp(theta)
              end
              ri = ri + x[j] * t
        end
        r[i] = yosb2[i] - ri
    end
    return r
end

function madsen_res(x, r)
    r[1] = x[1] ^ 2 + x[2] ^ 2 + x[1] * x[2]
    r[2] = sin(x[1])
    r[3] = cos(x[2])
    return r
end

function madsen_jac(x, j)
    j[1,1] = 2.e0 * x[1] + x[2]
    j[1,2] = 2.e0*x[2] + x[1]
    j[2,1] = cos(x[1])
    j[2,2] = 0.e0
    j[3,1] = 0.e0
    j[3,2] = -sin(x[2])
    return j
end

const ymeyer = [3.478e4, 2.861e4, 2.365e4, 1.963e4, 1.637e4, 1.372e4,
                1.154e4, 9.744e3, 8.261e3, 7.030e3, 6.005e3, 5.147e3,
                4.427e3, 3.820e3, 3.307e3, 2.872e3]

function meyer_res(x, r)
    for i = 1:16
        ti = Float64(5*i + 45)
        r[i] = x[1] * exp(x[2] / (ti + x[3])) - ymeyer[i]
    end
    return r
end

function meyer_jac(x, j)
    for i = 1:16
         ti = Float64(5*i + 45)
         u = ti + x[3]
         t = exp(x[2] / u)
         j[i, 1] = t
         j[i, 2] = x[1] * t/u
         j[i, 3] = -x[1] * x[2] * t / (u * u)
    end
    return j
end

# note for these problems m == p
function brown_res(x, r)
    n = length(x)
    t = x[1] - Float64(n + 1)
    for i = 2:n
        t = t + x[i]
    end
    nm1 = n - 1
    for i = 1:nm1
        r[i] = t + x[i]
    end
    t = x[1]
    for i = 2:n
        t = t * x[i]
    end
    r[n] = t - 1.0e0
    return r
end

brown5_res(x, r) = brown_res(x, r)
brown10_res(x, r) = brown_res(x, r)
brown30_res(x, r) = brown_res(x, r)
brown40_res(x, r) = brown_res(x, r)

function chebyquad_jac(x, j)
    n = size(j)[1]
    for k = 1:n
        tim1 = -1.e0 / Float64(n)
        z = 2.e0 * x[k] - 1.e0
        ti = z * tim1
        tpim1 = 0.e0
        tpi = 2.e0 * tim1
        z = z + z
        for i = 1:n
            j[i, k] = tpi
            tpip1 = 4.e0 * ti + z * tpi - tpim1
            tpim1 = tpi
            tpi = tpip1
            tip1 = z * ti - tim1
            tim1 = ti
            ti = tip1
        end
    end
    return j
end

function brown_dennis_jac(x, j)
    n = size(j)[1]
    for i = 1:n
         ti = 0.2e0 * Float64(i)
         j[i, 1] = 2.0e0 * (x[1] + x[2]*ti - exp(ti))
         j[i, 2] = ti * j[i, 1]
         t = sin(ti)
         j[i, 3] = 2.0e0 *(x[3] + x[4]*t - cos(ti))
         j[i, 4] = t * j[i, 3]
    end
    return j
end

function osborne2_jac(x, j)
    for i = 1:65
        ti = Float64(1 - i) * 1.e-1
        j[i, 1] = -exp(x[5] * ti)
        j[i, 5] = x[1] * ti * j[i,1]
        for k = 2:4
            t = x[k + 7] + ti
            theta = -x[k+4]*t*t
            theta > uftolg ? r2 = -exp(theta) : r2 = 0.e0
            j[i, k] = r2
            r2 = -t * r2 * x[k]
            j[i, k+4] = r2*t
            j[i, k+7] = 2.e0 * x[k+4] * r2
        end
    end
    return j
end

# note for these problems m == p
function brown_jac(x, j, n)
    nm1 = n - 1
    for k = 1:n
        for i = 1:nm1
            j[i, k] = 1.0e0
            if i == k
                j[i, k] = 2.0e0
            end
        end
    end
    for k = 1:n
        t = 1.0e0
        for i = 1:n
            i != k ? t = t * x[i] : t = t
        end
        j[n, k] = t
    end
    return j
end

brown5_jac(x, j) = brown_jac(x, j, 5)
brown10_jac(x, j) = brown_jac(x, j, 10)
brown30_jac(x, j) = brown_jac(x, j, 30)
brown40_jac(x, j) = brown_jac(x, j, 40)

nl2rj = Dict()

# Collect the residual and jacobian functions for nl2sol
for (k, v) in problems
    nl2rj[k] = (eval(Symbol(string(k, "_res"))), eval(Symbol(string(k, "_jac"))))
end


const cmplt = ["rosenbrock", "helix", "singular", "woods", "zangwill", 
               "engvall", "branin", "beale", "box", "freudenstein_roth",
               "watson6", "watson9", "watson12", "watson20", "chebyquad",
               "brown_dennis", "bard", "jennrich_sampson", "kowalik_osborne",
               "osborne1", "osborne2", "madsen", "meyer", "brown5", "brown10",
                "brown40"]

# Fortran test driver does not do [1, 10, 100] for these problems :
#   jennrich_sampson  [1.]   100*x_init will overflow R
#   osborne1   [1., 10.]
#   meyer  [1.] and must raise maxIter and maxFuncCall limits.  Get the same
#          answer as fortran driver but takes 192 iterations instead of 205
#   brown5, brown10,

const problematic =  ["cragg_levy", "davidon"]

#const working = cmplt[1:3]
const working = cmplt

# Set levenberg_marquardt x convergence tolerance in line with nl2sol. 
# Note that nl2sol does not directly test the convergence of the gradient
# as does levenberg_marquardt.  Rather it has tests for the relative and 
# absolute function convergence, so we don't mess with LM's tolG.
const tolX = sqrt(0.999*eps())
quiet = false

conv = Dict(
          3 => "x",
          4 => "r",
          5 => "b",
          6 => "a",
          7 => "s",
          8 => "f",
          9 => "e",
         10 => "i"
       )

function dump(res, fname)
    ios = open(fname, "w")
    fs = "{1:s}, {2:d}, {3:d}, {4:d}, {5:d}, {6:d}, {7:s}, {8:.2e}, {9:.2e}, {10:.2e}\n"
    sr = sortrows(res)
    for i in 1:size(sr)[1]
        write(ios, format(fs, replace(sr[i,1], r"_\d", ""), sr[i,2:end]...))
    end
    close(ios)
end

# To help in debugging NL2sol.  Only call on a single problem
function runone(;prb="rosenbrock", scale=1.0, mxfc=200, mxiter=200, sno=false)
    nlres, nljac = nl2rj[prb]
    r_init, j_init, x_init, s = problems[prb]
    n = length(r_init)
    p = length(x_init)
    iv, v = nl2_set_defaults(n, p)
    iv[MXFCAL] = mxfc
    iv[MXITER] = mxiter
    println("Running NL2sol on problem ", prb)
    results = nl2sol(nlres, nljac, scale * x_init, n, iv, v)
    println(return_code[iv[1]])
    println(results)

    if sno
        nl2_reset_defaults!(iv, v)
        iv[MXFCAL] = mxfc
        iv[MXITER] = mxiter
        println("Running NL2sno on problem ", prb)
        sno_results = nl2sno(nlres, scale * x_init, n, iv, v)
        println(return_code[iv[1]])
        println(sno_results)
    end

end

function runall()
    # So it turns out that if we run the nl2sno tests, then we will crash
    # deep, deep within garbage collecion or get a corrupted double linked list
    # message from glibc.  If we disable garbage collection and use the newly
    # added nl2_reset_defaults! for nl2sno instead of allocating new ones, then
    # this set of tests will run to completion about 1/2 the time.
    use_nl2sno && gc_enable(false)
    all_results = nothing
    for (prb, v) in problems
        nlres, nljac = nl2rj[prb]
        lmres, lmjac = lmrj[prb]
        r, j, x, s = v
        n = length(r)
        p = length(x)
        scale = 1.0
        for i in 1:s
            x_init = scale * x
            nl_results = "No NL2sol results available"
            results = "No LM results available"
            sno_results = "No NL2sno results available"
            iv, v = nl2_set_defaults(n, p)
            # meyer needs larger limits
            if prb == "meyer"
                iv[MXFCAL] = 350
                iv[MXITER] = 350
            elseif prb == "watson20"
                # match the fortran test driver
                iv[MXITER] = 15
            end
            iv[PRUNIT] = 0  # supress nl2sol output
            
            nl_results = try
                println("\nStarting NL2sol on problem $prb at scale $scale")
                nl2sol(nlres, nljac, x_init, n, iv, v)
            catch exc
                println(exc)
                println("NL2sol exception $exc on problem $prb")            
            end

            # # We subtract 1 from i (the scaling) to match the NL2sol paper
            tag = string(prb, "_", i-1)

            # t is the predicted relative function reduction.  Note that
            # in Table II of the paper it is always positive.  However, 
            # we find a few negative values.  These match the full Fortran 
            # version when that is run, so I suspect that they only printed
            # the abs() of this value in Table II, so that is what I do here
            v[FUNCT0] > 0.0 ? t = abs(v[NREDUC] / v[FUNCT0]) : t = 1.0

            # This summary line is meant to match "Table II Default NL2SOL"
            # in "An Adaptive Nonlinear Least-Squares Algorithm", ie the
            # algorithm paper.
            # Remember to subtract out the number of function calls
            # and gradient calls used to calculate the covarience..
            #
            # NOTE: some new (Version 0.5.0-dev+5332) weird issue causes
            # problems with arrays of type string so that I cannot make a vector
            # of strings and then take the transpose.  To make a row, I need a
            # very long line... check later to see if this gets fixed.
            nl = [tag (i-1)  n  p  Int(iv[NFCALL] - iv[NFCOV])  Int(iv[NGCALL] - iv[NGCOV])  conv[Int(iv[1])] v[FUNCT]  t   v[RELDX]]

            if all_results == nothing
                all_results = nl
            else
                all_results = [all_results; nl]
            end
            
            results = try
                if use_levenberg
                    println("\nStarting Levenberg Marquardt on problem  $prb at scale $scale")
                    levenberg_marquardt(lmres, lmjac, x_init; 
                                        maxIter=400, tolX=tolX)
                end
            catch exc
                println(exc)
                println("Levenberg Marquardt exception $exc on problem $prb")
            end

            if use_nl2sno
                nl2_reset_defaults!(iv, v)
                iv[PRUNIT] = 0
                sno_results = try
                    println("\nStarting NL2sno on problem  $prb at scale $scale")
                    nl2sno(nlres, x_init, n, iv, v)
                catch exc
                    println(exc)
                    println("NL2sno exception $exc on problem $prb")
                end
            end

            if !quiet
                println("\nnl2sol on problem $prb at scale $scale")
                println(nl_results)
                if use_levenberg
                    println("\nlevenberg-marquardt on problem $prb at scale $scale")
                    println(convert_results(results))
                end
                if use_nl2sno
                    println("\nNL2sno on problem $prb at scale $scale")
                    println(sno_results)
                end
            end

            scale *= 10.0
        end
    end
    dump(all_results, "nlresults.log")

    ! use_DataFrames && return true
    
    origDF = readtable("nl2_results.txt", header=false)
    newDF  = readtable("nlresults.log", header=false)
    pass = true

    # compare just a few results now.  Not all the parameters, or
    # even the problem names, are the same between the Julia versions
    # and the Fortran versions.  It goes without saying that this is
    # incredibly brittle.
    check = [(1,1), (5,5), (25, 29), (26, 30), (42, 46), (43, 47), (44, 48)]
    for (i, j) in check
        origDF[i, :] != newDF[j, :] ? pass = false : nothing
    end
    println("Passed subset of NL2sol paper results")
    #gc_enable(true)  # we get farther using nl2sno, generally, if we don't re-enable.  Why?
    return pass
end

# There is a lot of noise from levenberg_marquardt, but oh well for now
!isinteractive() && @test(runall())
