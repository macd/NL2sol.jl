module NL2sol
import Base
using Printf
using Random
using LinearAlgebra
using Pkg

export nl2sol, nl2sno, nl2_set_defaults, nl2_reset_defaults!, return_code
export MXFCAL, MXITER, OUTLEV, PRUNIT, NFCALL, NGCALL, NITER, NFCOV, NGCOV
export AFCTOL, RFCTOL, XCTOL, XFTOL, NREDUC, DGNORM, DSTNRM, PREDUC, RADIUS
export FUNCT, FUNCT0, RELDX

# OK, call me paranoid, but lets add some padding for NL2SOL
const PADDING = 1000

# Borrowed from Optim so we don't have a dependency
mutable struct NL2OptimizationResults{T, N}
    method::String
    initial_x::Array{T, N}
    minimum::Array{T, N}
    f_minimum::Float64
    iterations::Int
    iteration_converged::Bool
    x_converged::Bool
    x_tol::Float64
    f_converged::Bool
    f_tol::Float64
    g_converged::Bool
    g_tol::Float64
    f_calls::Int
    g_calls::Int
end

function Base.show(io::IO, r::NL2OptimizationResults)
    @printf io "Results of Optimization Algorithm\n"
    @printf io " * Algorithm: %s\n" r.method

    if length(join(r.initial_x, ",")) < 40
        @printf io " * Starting Point: [%s]\n" join(r.initial_x, ",")
    else
        @printf io " * Starting Point: [%s, ...]\n" join(r.initial_x[1:2], ",")
    end
    if length(join(r.minimum, ",")) < 40
        @printf io " * Minimizer: [%s]\n" join(r.minimum, ",")
    else
        @printf io " * Minimizer: [%s, ...]\n" join(r.minimum[1:2], ",")
    end
    @printf io " * Minimum: %e\n" r.f_minimum
    @printf io " * Iterations: %d\n" r.iterations
    @printf io " * Convergence: %s\n" r.x_converged || r.f_converged || r.g_converged
    @printf io "   * |x - x'| < %.1e: %s\n" r.x_tol r.x_converged
    @printf io "   * |f(x) - f(x')| / |f(x)| < %.1e: %s\n" r.f_tol r.f_converged
    @printf io "   * |g(x)| < %.1e: %s\n" r.g_tol r.g_converged
    @printf io "   * Reached Maximum Number of Iterations: %s\n" r.iteration_converged
    @printf io " * Objective Function Calls: %d\n" r.f_calls
    @printf io " * Gradient Calls: %d\n" r.g_calls
    return
end

# NL2SOL Return codes.  Returned in iv[1]
const return_code = 
    Dict(3  => "x convergence",
         4  => "relative function convergence",
         5  => "both x and relative function convergence",
         6  => "absolute function convergence",
         7  => "singular convergence",
         8  => "false convergence",
         9  => "function evaluation limit",
         10 => "iteration limit",
         11 => "internal stopx",
         13 => "f(x) cannot be computed at the initial x",
         14 => "Bad parameters passed to assess",
         15 => "The Jacobian could not be computed at x",
         16 => "n or p out of range: p < 0 or n < p",
         17 => "A restart was attempted with n, p changed"
    )

# Input values stored in iv
const MXFCAL = 17
const MXITER = 18
const OUTLEV = 19
const PRUNIT = 21

# Output values stored in iv
const NFCALL =  6
const NGCALL = 30
const NITER  = 31
const NFCOV  = 40
const NGCOV  = 41

# Input values in v controlling convergence
const AFCTOL = 31   # absolute function convergence default = max(10e-20, eps^2)
const RFCTOL = 32   # relative function convergence default = max(10e-10, eps^(2/3))
const XCTOL  = 33   # x convergence default = sqrt(eps)
const XFTOL  = 34   # false convergence default = 100*eps

# Output values stored in v
const DGNORM =  1
const DSTNRM =  2
const NREDUC =  6
const PREDUC =  7
const RADIUS =  8
const FUNCT  = 10
const FUNCT0 = 13
const RELDX  = 17

# Default values for convergence testing
const df_tolX = sqrt(0.999 * eps())  
const df_tolFalseX = 100*eps()
const df_tolAbsFunc = max(1e-20, eps()^2)
const df_tolRelFunc = max(1e-10, eps()^(2/3))
const df_maxIter = 400
const df_maxFuncCall = 400

# Note that the lib path needs to be a constant string
if haskey(ENV, "NL2SOL_LIBPATH")
    const libnl2sol = joinpath(ENV["NL2SOL_LIBPATH"], "libnl2sol.so")
else
    const libnl2sol = joinpath(@__DIR__, "../deps/usr/lib/libnl2sol.so")
end

function nl2_reset_defaults!(iv, v)
    iv[:] = 0
    v[:] = 0.0
    ccall((:dfault_, libnl2sol), Cvoid, (Ref{Int32}, Ref{Float64}), iv, v)
end

function nl2_set_defaults(n, p)
    ivsize = p + 60 + PADDING
    vsize = ceil(Int, 93 + n*(p + 3) + (3 * p * (p + 11))/2) + PADDING
    iv = zeros(Int32, ivsize)
    v  = zeros(Float64, vsize)
    ccall((:dfault_, libnl2sol), Cvoid, (Ref{Int32}, Ref{Float64}), iv, v)
    return iv, v
end

function nl2_set_residual(res::Function)
    # Have tried function cacheing here in the past but probably too 
    # dangerous for by chance (or programmtically generated names) may get
    # same name with different definition.  Also, if we try to redefine even
    # with a new function definition,
    # commit 89424cc05a3fae94221efc45f24f924a75d2f58a
    # causes some kind of corruption on the call to the redefined function 
    # (not the first one), so just make certain to have a unique name now.
    wr = Symbol(string("nl2_", res, "_", randstring(5)))
    # Wrap residual for the nl2sol/nl2sno calling signature
    func = quote
        function ($wr)(n_, p_, x_::Ref{Float64}, nf_::Ref{Int32}, r_::Ref{Float64}, 
                          uiparm, urparm, ufparm)
            n = unsafe_load(n_, 1)
            p = unsafe_load(p_, 1)
            x = unsafe_wrap(Array, x_, p)
            r = unsafe_wrap(Array, r_, n)

            # If the residual calculation raises a DomainError, we have taken 
            # a step inside of NL2SOL that is too big.  By setting nf_ to zero,
            # we are telling it to take a smaller step
            try
                ($res)(x, r)
            catch y
                if isa(y, DomainError)
                    unsafe_store!(nf_, Int32(0))
                else
                    throw(y)
                end
            end
            return
        end
    end
    nlf = eval(func)

    # Using the new 0.7 / 1.0 form as a macro results in nlf not being
    # defined unless we interpolate it into the macro call.  This is
    # discussed in the manual in the chaper "Calling C and Fortran Code"
    # in section "Closure cfunctions".
    # I also needed to use Ptr{...} rather than Ref{...}, otherwise we get
    # a hard crash.
    nr = @cfunction($nlf, Cvoid, (Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, 
                                  Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, 
                                  Ptr{Float64}, Ptr{Ptr{Cvoid}}))

    return nr
end

function nl2_set_jacobian(jacobian::Function)
    # See comments above on names
    wj = Symbol(string("nl2_", jacobian, "_", randstring(5)))
        
    # Wrap jacobian for nl2sol calling signature
    nj = quote
        function ($wj)(n_::Ptr{Int32},
                       p_::Ptr{Int32},
                       x_::Ref{T},
                       nf_::Ptr{Int32},
                       jac_::Ref{T},
                       uiparm::Ptr{Int32},
                       urparm::Ptr{Float64},
                       ufparm::Ptr{Ptr{Nothing}}) where {T}
            n = unsafe_load(n_, 1)
            p = unsafe_load(p_, 1)
            x = unsafe_wrap(Array, x_, p)
            jac = unsafe_wrap(Array, jac_, (n, p))
            ($jacobian)(x, jac)
            return
        end
    end
    nlj = eval(nj)

    # see comments above
    jc = @cfunction($nlj, Nothing, (Ptr{Int32}, Ptr{Int32}, Ptr{Float64},
                                    Ptr{Int32}, Ptr{Float64}, Ptr{Int32},
                                    Ptr{Float64}, Ptr{Ptr{Nothing}}))
    
    return jc
end

mutable struct NL2Results{T}
    algorithm::AbstractString
    initial_x::T
    final_x::T
    rcode::AbstractString
end

function nl2sno(res::Function, init_x, n, iv, v)
    p = length(init_x)
    x = copy(init_x)

    # These need to be 32 bit ints for nl2sol
    p_ = Int32(p)   
    n_ = Int32(n)

    # Currently, we do not use any of the u*parm arrays.  (They are only
    # in NL2SOL because Fortran did not have closures)
    uiparm = Int32[]
    urparm = Float64[]
    ufparm = Array{Ref{Cvoid}}(undef, 1)
    nl2res = nl2_set_residual(res)

    ccall((:nl2sno_, libnl2sol), Nothing,
        (Ref{Int32},
         Ref{Int32},
         Ref{Float64},
         Ptr{Nothing},  # this MUST be a Ptr{Nothing}
         Ref{Int32},
         Ref{Float64},
         Ref{Int32},
         Ref{Float64},
         Ref{Ref{Cvoid}}),
         n_, p_, x, nl2res, iv, v, uiparm, urparm, ufparm)

    (iv[end] != 0 || v[end] != 0.0) && error("NL2SNO memory corruption")

    results = NL2OptimizationResults(
        "nl2sno",
         init_x,
         x,
         v[FUNCT],
         Int(iv[NITER]),
         Int(iv[NITER]) >= Int(iv[MXFCAL]),
         Int(iv[1]) == 3 || Int(iv[1]) == 5,
         v[RELDX],
         Int(iv[1]) == 4 || Int(iv[1]) == 6,
         0.0, 
         false,
         0.0,
         # TODO: check these against paper
         Int(iv[NFCALL] - iv[NFCOV]),
         Int(iv[NGCALL] - iv[NGCOV])
    )
    return results
end


"""
    nl2sol solves the non-linear least squares problem.  That is, it finds
    an x that minimizes  sum(i=1:n){r_i(x)^2} where x is a vector of size p.
    It returns a struct of type Optim.MultivariateOptimizationResults that
    contains the relevant info (see the Optim docs for further info)

    The residual and the jacobian functions are expected to take args that
    have been preallocated for those values.  These arrays are actually 
    allocated in the Julia function nl2sol before passing to the Fortran
    subroutine nl2sol.

    NOTE: NL2sol.jl does pointer tricks in order to interface with the
          FORTRAN code, so tread lightly if you'd like to modify the code.

    EXAMPLE USAGE

    using NL2sol

    function rosenbrock_res(x, r)
        r[1] = 10. * (x[2] - x[1]^2 )
        r[2] = 1. - x[1]
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

"""
function nl2sol(res::Function, jac::Function, init_x, n, iv, v)
    p = length(init_x)
    x = copy(init_x)
    p_ = Int32(p)
    n_ = Int32(n)

    # Currently, we do not use any of the u*parm arrays.
    uiparm = Array{Int32}(undef, 1)
    urparm = Float64[]
    ufparm = Array{Ref{Cvoid}}(undef, 1)

    nl2res = nl2_set_residual(res)
    nl2jac = nl2_set_jacobian(jac)

    ccall((:nl2sol_, libnl2sol), Cvoid,
        (Ref{Int32},
         Ref{Int32},
         Ref{Float64},
         Ptr{Nothing},   # this MUST be a Ptr{Nothing}
         Ptr{Nothing},   # this MUST be a Ptr{Nothing}
         Ref{Int32},
         Ref{Float64},
         Ref{Int32},
         Ref{Float64},
         Ref{Ref{Cvoid}}),
         n_, p_, x, nl2res, nl2jac, iv, v, uiparm, urparm, ufparm)

    (iv[end] != 0 || v[end] != 0.0) && error("NL2SOL memory corruption")

    results = NL2OptimizationResults(
        "nl2sol",
         init_x,
         x,
         v[FUNCT],
         Int(iv[NITER]),
         Int(iv[NITER]) >= Int(iv[MXFCAL]),
         Int(iv[1]) == 3 || Int(iv[1]) == 5,
         v[RELDX],
         Int(iv[1]) == 4 || Int(iv[1]) == 6,
         0.0, 
         false,
         0.0,
         Int(iv[NFCALL] - iv[NFCOV]),
         Int(iv[NGCALL] - iv[NGCOV])
                                  )
    return results
end

# Convenience function for nl2sol.  Don't have to mess with iv and v but
# can still override some of the most common ways to control an optimization
# algorithm.
function nl2sol(res::Function, jac::Function, init_x, n; 
                maxIter=df_maxIter, maxFuncCall=df_maxFuncCall, 
                tolX=df_tolX,  tolAbsFunc=df_tolAbsFunc,
                tolRelFunc=df_tolRelFunc, quiet=true)
    p = length(init_x)
    iv, v = nl2_set_defaults(n, p)
    # Set defaults
    iv[MXITER] = maxIter
    iv[MXFCAL] = maxFuncCall
    v[XCTOL] = tolX
    v[AFCTOL] = tolAbsFunc
    v[RFCTOL] = tolRelFunc
    quiet ? iv[PRUNIT] = 0 : nothing
    results = nl2sol(res, jac, init_x, n, iv, v)
    quiet || println("Convergence state: ", return_code[iv[1]])
    return results
end

end # module NL2sol
