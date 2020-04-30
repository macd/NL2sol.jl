const zangwill = let res_init=zeros(3), jac_init=zeros(3,3), x_init=[100., -1., 2.5]

    function res(x, r)
        r[1] = x[1] - x[2] + x[3]
        r[2] = -x[1] + x[2] + x[3]
        r[3] = x[1] + x[2] - x[3]
        return r
    end

    function jac(x, jac)
        jac[:] .= 0.0
        for k = 1:3, i = 1:3
            jac[i, k] = 1.0
        end
        jac[1, 2] = -1.0
        jac[2, 1] = -1.0
        jac[3, 3] = -1.0
        return jac
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("zangwill", res, jac,
                                                   res_init, jac_init, x_init;
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end
