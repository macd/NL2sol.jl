const engvall = let res_init=zeros(5), jac_init=zeros(5,3), x_init=[1.0, 2.0, 0.0]

    function res(x, r)
        r[1] = x[1]^2 + x[2]^2 + x[3]^2 - 1.0
        r[2] = x[1]^2 + x[2]^2 + (x[3] - 2.0)^2 - 1.0
        r[3] = x[1] + x[2] + x[3] - 1.0
        r[4] = x[1] + x[2] - x[3] + 1.0
        r[5] = x[1]^3 + 3.0 * x[2] ^ 2 + (5.0 * x[3] - x[1] + 1.0)^2 - 36.0
        return r
    end

    function jac(x, jac)
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

    f(;scale=1, verbose=false, print_steps=false) = testone("engvall", res, jac,
                                                   res_init, jac_init, x_init;
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end
