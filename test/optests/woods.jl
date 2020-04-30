const woods = let res_init=zeros(7), jac_init=zeros(7,4), x_init=[-3.0, -1.0, -3.0, -1.0]

    function res(x, r)
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

    function jac(x, jac)
        jac[:] .= 0.0
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

    f(;scale=1, verbose=false, print_steps=false) = testone("woods", res, jac,
                                                   res_init, jac_init, x_init;
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end
