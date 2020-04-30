const freudenstein_roth = let res_init=zeros(2), jac_init=zeros(2,2), x_init=[15.0, -2.0]

    function res(x, r)
        r[1] = -1.3e1 + x[1] - 2.0e0*x[2] + 5.0e0*x[2]^2 - x[2]^3
        r[2] = -2.9e1 + x[1] - 1.4e1*x[2] + x[2]^2 + x[2]^3
        return r
    end

    function jac(x, j)
        j[1,1] = 1.e0
        j[1,2] = -2.e0 + x[2]*(1.e1 - 3.e0*x[2])
        j[2,1] = 1.e0
        j[2,2] = -1.4e1 + x[2]*(2.e0 + 3.e0*x[2])
        return j
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("freudenstein_roth", res, jac,
                                                   res_init, jac_init, x_init; 
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end
