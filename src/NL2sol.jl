module NL2sol
import Base
import Optim

using Lexicon
using Docile
@document

export nl2sol

# OK, call me paranoid, but lets add some padding for NL2SOL
const PADDING = 1000

# NL2SOL Return codes
const return_code = 
    Dict{Int, ASCIIString}(3  => "x convergence",
         4  => "relative function convergence",
         5  => "both x and relative function convergence",
         6  => "absolute function convergence",
         7  => "singular convergence",
         8  => "false convergence",
         9  => "function evaluation limit",
         10 => "iteration limit",
         11 => "internal stopx"
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
const RADIUS =  8
const FUNCT  = 10
const RELDX  = 17

# Default values for convergence testing
const df_tolX = sqrt(0.999 * eps())  
const df_tolFalseX = 100*eps()
const df_tolAbsFunc = max(1e-20, eps()^2)
const df_tolRelFunc = max(1e-10, eps()^(2/3))
const df_maxIter = 400
const df_maxFuncCall = 400

# Note that the lib needs to be a constant
if haskey(ENV, "NL2SOL_LIBPATH")
    const libnl2sol = joinpath(ENV["NL2SOL_LIBPATH"], "libnl2sol.so")
else
    const libnl2sol = joinpath(Pkg.dir(), "NL2sol/deps/usr/lib/libnl2sol.so")
end

function nl2sol_set_defaults(iv, v)
    ccall((:dfault_, libnl2sol), Void, (Ptr{Int32}, Ptr{Float64}), iv, v)
end

type NL2Array{T}
    p::Ptr{T}
    rows::Int32
end

Base.getindex(x::NL2Array, i) = unsafe_load(x.p, i)
Base.setindex!(x::NL2Array, y, i) = unsafe_store!(x.p, y, i)
Base.length(x::NL2Array) = x.rows
Base.endof(x::NL2Array) = length(x)
Base.start(::NL2Array) = 1    #start(x::NL2Array) = 1
Base.next(x::NL2Array, i) = (x[i], i+1)
Base.done(x::NL2Array, i) = (i > length(x))

type NL2Matrix{T}
    p::Ptr{T}
    rows::Int32
    cols::Int32
end

Base.getindex(x::NL2Matrix, i, j) = unsafe_load(x.p, x.rows*(j - 1) + i)
Base.getindex(x::NL2Matrix, i) = unsafe_load(x.p, i)  # as one-D
Base.setindex!(x::NL2Matrix, y, i, j) = unsafe_store!(x.p, y, x.rows*(j - 1) + i)
Base.setindex!(x::NL2Matrix, y, i) = unsafe_store!(x.p, y, i) # as one-D

# living on the edge...  This is to make x[:] = 0.0 work for the NL2 types
Base.unsafe_store!(x::Ptr{Float64}, val::Float64, ::Colon) =
    function unsafe_store!(x::Ptr{Float64}, val::Float64, ::Colon)
        l = length(x)
        for i = 1:l
            unsafe_store!(x, val, i)
        end
    end

Base.unsafe_store!(x::Ptr{Float64}, val::Float64, I::UnitRange{Int}) =
    function unsafe_store!(x::Ptr{Float64}, val::Float64, I::UnitRange{Int})
        for i in I
            unsafe_store!(x, val, i)
        end
    end

Base.length(x::NL2Matrix) = x.rows * x.cols
Base.endof(x::NL2Matrix) = length(x)
Base.size(x::NL2Matrix) = (x.rows, x.cols)

function nl2sol_set_functions(res, jac)
    wr = Symbol(string("nl2_", res))
    cr = Symbol(string(wr, "_cr"))
    wj = Symbol(string("nl2_", jac))
    cj = Symbol(string(wj, "_cj"))
    @eval begin
        # Wrap residual for nl2sol calling signature
        function ($wr){T}(n_, p_, x_::Ptr{T}, nf_, r_::Ptr{T}, uiparm, urparm, ufparm)
            n = unsafe_load(n_, 1)
            p = unsafe_load(p_, 1)
            x = NL2Array(x_, p)
            r = NL2Array(r_, n)
            ($res)(x, r)
            return
        end
        # Now make it C (actually Fortran) callable
        const $(cr) = cfunction(($wr), Void, 
                                (Ptr{Int32},
                                 Ptr{Int32}, 
                                 Ptr{Float64}, 
                                 Ptr{Int32}, 
                                 Ptr{Float64}, 
                                 Ptr{Int32}, 
                                 Ptr{Float64}, 
                                 Ptr{Ptr{Void}}))

        # Wrap jacobian for nl2sol calling signature
        function ($wj){T}(n_, p_, x_::Ptr{T}, nf_, jac_::Ptr{T}, uiparm, urparm, ufparm)
            n = unsafe_load(n_, 1)
            p = unsafe_load(p_, 1)
            x = NL2Array(x_, p)
            jac = NL2Matrix(jac_, n, p)
            ($jac)(x, jac)
            return
        end
        # and make it C callable
        const $(cj) = cfunction($(wj), Void, 
                                (Ptr{Int32}, 
                                 Ptr{Int32}, 
                                 Ptr{Float64}, 
                                 Ptr{Int32}, 
                                 Ptr{Float64}, 
                                 Ptr{Int32}, 
                                 Ptr{Float64}, 
                                 Ptr{Ptr{Void}}))
    end
    return (eval(cr), eval(cj))
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

    NOTE: NL2sol does nasty pointer tricks that are liable to make programs
          brittle, so don't do a "using NL2sol" in code you care about.

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
function nl2sol(res::Function, jac::Function, init_x, n; 
                tolX=df_tolX,  tolFalseX=df_tolFalseX, tolAbsFunc=df_tolAbsFunc,
                tolRelFunc=df_tolRelFunc, maxIter=df_maxIter, 
                maxFuncCall=df_maxFuncCall, quiet=true)
    p = length(init_x)
    x = copy(init_x)
    p_ = Int32(p)
    n_ = Int32(n)

    ivsize = p + 60 + PADDING
    vsize = round(Int, 93 + n*(p + 3) + (3 * p * (p + 11))/2) + PADDING
    iv = zeros(Int32, ivsize)
    v  = zeros(Float64, vsize)
    nl2sol_set_defaults(iv, v)
    # Set defaults
    iv[MXITER] = maxIter
    iv[MXFCAL] = maxFuncCall
    v[XCTOL] = tolX
    v[XFTOL] = tolFalseX
    v[AFCTOL] = tolAbsFunc
    v[RFCTOL] = tolRelFunc
    quiet ? iv[PRUNIT] = 0 : nothing
    # Currently, we do not use any of these.
    uiparm = Array(Int32, 1)
    urparm = Float64[]
    ufparm = Array(Ptr{Void}, 1)
    nl2res, nl2jac = nl2sol_set_functions(res, jac)

    ccall((:nl2sol_, libnl2sol), Void,
        (Ptr{Int32},    # many need to change this to Tuple{Type1, Type2,...}
         Ptr{Int32},
         Ptr{Float64},
         Ptr{Void},
         Ptr{Void},
         Ptr{Int32},
         Ptr{Float64},
         Ptr{Int32},
         Ptr{Float64},
         Ptr{Void}),
         &n_, &p_, x, nl2res, nl2jac, iv, v, uiparm, urparm, ufparm)

    (iv[end] != 0 || v[end] != 0.0) && error("NL2SOL memory corruption")

    quiet || println("Convergence state: ", return_code[iv[1]])
    results = Optim.MultivariateOptimizationResults(
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
         Optim.OptimizationTrace(),
         Int(iv[NFCALL] - iv[NFCOV]),
         Int(iv[NGCALL] - iv[NGCOV])
    )
    return results
end

end # module NL2sol
