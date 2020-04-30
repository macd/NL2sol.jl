const branin = let res_init=zeros(2), jac_init=zeros(2,2), x_init=[2.0, 0.0]

    function res(x, r)
        r[1] = 4.0 * (x[1] + x[2])
        r[2] = r[1] + (x[1] - x[2]) * ((x[1] - 2.0) ^ 2 + x[2] ^ 2 - 1.0)
        return r
    end

    function jac(x, jac)
        jac[1, 1] = 4.0
        jac[1, 2] = 4.0
        jac[2, 1] = 3.0 + (x[1] - 2.0) * (3.0 * x[1] - 2.0 * x[2] - 2.0) + x[2] * x[2]
        jac[2, 2] = 1.0 + 2.0 * (2.0*x[1] - x[2] * x[2]) - (x[1] - x[2])^2
        return jac
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("branin", res, jac,
                                                            res_init, jac_init, x_init;
                                                            scale=scale, verbose=verbose,
                                                            print_steps=print_steps)
end
