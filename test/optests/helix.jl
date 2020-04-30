const helix = let res_init=zeros(3), jac_init=zeros(3,3), x_init=[-1.0, 0.0, 0.0]

    function res(x, r)
        theta = atan(x[2], x[1]) / 2pi

        if x[1] <= 0.0 && x[2] <= 0.0
            theta = theta + 1.0
        end

        r[1] = 10.0 * (x[3] - 10.0 * theta )
        r[2] = 10.0 * (sqrt(x[1]^2 + x[2]^2) - 1.0)
        r[3] = x[3]
        return r
    end

    function jac(x, jac)
        t = x[1]^2 + x[2]^2
        ti = 100.0 / (2pi * t )
        jac[1, 1] = ti * x[2]
        t = 10.0 / sqrt(t)
        jac[2, 1] = x[1] * t
        jac[3, 1] = 0.0
        jac[1, 2] = -ti * x[1]
        jac[2, 2] = x[2] * t
        jac[3, 2] = 0.0
        jac[1, 3] = 10.0
        jac[2, 3] = 0.0
        jac[3, 3] = 1.0
        return jac
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("helix", res, jac,
                                                            res_init, jac_init, x_init;
                                                            scale=scale, verbose=verbose,
                                                            print_steps=print_steps)
end
