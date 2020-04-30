# It turns out that this works OK for starting points 10*x_init and
# 100*x_init.  Interesting, the Fortran code does _not_ run the case
# 1*x_init.   It looks like it converges to a different minimum from
# that starting point.
const expmax = 1.999 * log(sqrt(0.999*floatmax()))
const expmin = 1.999 * log(sqrt(1.001*floatmin()))

const box = let res_init=zeros(10), jac_init=zeros(10,3), x_init=[0.0, 1.0, 2.0]

    function res(x, r)
        -expmax > minimum(x) && throw("Error in box residual")
        for i = 1:10
            ti = -0.1 * Float64(i)
            t1 = ti * x[1]
            e1 = t1 > expmin ? exp(t1) : 0.0
            t2 = ti*x[2]
            e2 = t2 > expmin ? exp(t2) : 0.0
            r[i] = (e1 - e2) - x[3]*(exp(ti) - exp(10.0 * ti))
        end
        return r
    end

    function jac(x, jac)
        for i = 1:10
            ti = -0.1 * Float64(i)
            t = x[1] * ti
            e = t > expmin ? exp(t) : 0.0
            jac[i, 1] = ti * e
            t = x[2] * ti
            e = t > expmin ? exp(t) : 0.0
            jac[i, 2] = -ti * e
            jac[i, 3] = exp(10.0 * ti) - exp(ti)
        end
        return jac
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("box", res, jac,
                                                   res_init, jac_init, x_init;
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end
