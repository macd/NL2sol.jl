const beale = let res_init=zeros(3), jac_init=zeros(3,2), x_init=[0.1, 0.1]

    function res(x, r)
        r[1] = 1.5 - x[1] * (1.0 - x[2])
        r[2] = 2.25 - x[1] * (1.0 - x[2]^2)
        r[3] = 2.625 - x[1] * (1.0 -  x[2]^3)
        return r
    end

    function jac(x, jac)
        jac[1, 1] = x[2] - 1.0
        jac[1, 2] = x[1]
        jac[2, 1] = x[2] ^ 2 - 1.0
        jac[2, 2] = 2.0 * x[1] * x[2]
        jac[3, 1] = x[2] ^ 3 - 1.0
        jac[3, 2] = 3.0 * x[1] * (x[2] ^ 2)
        return jac
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("beale", res, jac,
                                                   res_init, jac_init, x_init;
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end
