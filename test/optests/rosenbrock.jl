const rosenbrock = let res_init=zeros(2), jac_init=zeros(2,2), x_init=[-1.2, 1.0]

    function res(x, r)
        r[1] = 10.0 * (x[2] - x[1]^2 )
        r[2] = 1.0 - x[1]
        return r
    end

    function jac(x, jac)
        jac[1, 1] = -20.0 * x[1]
        jac[1, 2] =  10.0
        jac[2, 1] =  -1.0
        jac[2, 2] =   0.0
        return jac
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("rosenbrock", res, jac,
                                                            res_init, jac_init, x_init;
                                                            scale=scale, verbose=verbose,
                                                            print_steps=print_steps)
end
