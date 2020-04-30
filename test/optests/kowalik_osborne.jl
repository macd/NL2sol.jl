const ykow = [1.957e-1, 1.947e-1, 1.735e-1, 1.600e-1, 8.44e-2, 6.27e-2,
              4.56e-2,  3.42e-2, 3.23e-2, 2.35e-2, 2.46e-2]

const ukow = [4.0e0, 2.0e0, 1.0e0, 5.0e-1, 2.5e-1, 1.67e-1,
              1.25e-1, 1.0e-1, 8.33e-2, 7.14e-2, 6.25e-2]

const kowalik_osborne = let res_init=zeros(11), jac_init=zeros(11,4), x_init=[.25, .39, .415, .39]

    function res(x, r)
        for i = 1:11
            r[i] = ykow[i] - x[1]*(ukow[i] ^ 2 + x[2] * ukow[i]) / (ukow[i] ^ 2 + x[3]*ukow[i] + x[4])
        end
        return r
    end

    function jac(x, j)
        for i = 1:11
            t = -1.e0 / (ukow[i]^2 + x[3] * ukow[i] + x[4])
            j[i, 1] = t * (ukow[i]^2 + x[2] * ukow[i])
            j[i, 2] = x[1] * ukow[i]*t
            t = t * j[i, 1] * x[1]
            j[i, 3] = ukow[i] * t
            j[i, 4] = t
        end
        return j
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("kowalik_osborne", res, jac,
                                                            res_init, jac_init, x_init;
                                                            scale=scale, verbose=verbose,
                                                            print_steps=print_steps)
end
