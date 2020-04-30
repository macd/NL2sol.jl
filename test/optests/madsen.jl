const madsen = let res_init=zeros(3), jac_init=zeros(3,2), x_init=[3.0, 1.0]

    function res(x, r)
        r[1] = x[1] ^ 2 + x[2] ^ 2 + x[1] * x[2]
        r[2] = sin(x[1])
        r[3] = cos(x[2])
        return r
    end

    function jac(x, j)
        j[1,1] = 2.e0 * x[1] + x[2]
        j[1,2] = 2.e0*x[2] + x[1]
        j[2,1] = cos(x[1])
        j[2,2] = 0.e0
        j[3,1] = 0.e0
        j[3,2] = -sin(x[2])
        return j
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("madsen", res, jac,
                                                   res_init, jac_init, x_init;
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end
