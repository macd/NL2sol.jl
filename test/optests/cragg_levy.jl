const cragg_levy = let res_init=zeros(5), jac_init=zeros(5,4), x_init=[1.0,2.0,2.0,2.0]

    function res(x, r)
        r[1] = (exp(x[1]) - x[2])^2
        r[2] = 10.0 * (x[2] - x[3])^3
        r[3] = (sin(x[3] - x[4]) / cos(x[3] - x[4]) )^2
        r[4] = x[1]^4
        r[5] = x[4] - 1.0
        return r
    end

    function jac(x, jac)
        jac[:] .= 0.0
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

    f(;scale=1, verbose=false, print_steps=false) = testone("cragg_levy", res, jac,
                                                   res_init, jac_init, x_init;
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end
