const singular = let res_init=zeros(4), jac_init=zeros(4,4), x_init=[3.0, -1.0, 0.0, 1.0]

    function res(x, r)
        r[1] = x[1] + 10.0 * x[2]
        r[2] = sqrt(5.0) * (x[3] - x[4])
        r[3] = (x[2] - 2x[3])^2
        r[4] = sqrt(10.0) * (x[1] - x[4])^2
        return r
    end

    function jac(x, jac)
        jac[:] .= 0.0
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

    f(;scale=1, verbose=false, print_steps=false) = testone("singular", res, jac,
                                                   res_init, jac_init, x_init;
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end

