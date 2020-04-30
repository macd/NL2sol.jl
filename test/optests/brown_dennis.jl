const brown_dennis = let res_init=zeros(20), jac_init=zeros(20,4), x_init=[25., 5., -5., -1.]

    function res(x, r)
        n = length(r)
        for i = 1:n
            ti = 0.2e0 * Float64(i)
            r[i] = (x[1] + x[2] * ti - exp(ti)) ^ 2 +
                (x[3] + x[4] * sin(ti) - cos(ti))^2
        end
        return r
    end

    function jac(x, j)
        n = size(j)[1]
        for i = 1:n
            ti = 0.2e0 * Float64(i)
            j[i, 1] = 2.0e0 * (x[1] + x[2]*ti - exp(ti))
            j[i, 2] = ti * j[i, 1]
            t = sin(ti)
            j[i, 3] = 2.0e0 *(x[3] + x[4]*t - cos(ti))
            j[i, 4] = t * j[i, 3]
        end
        return j
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("brown_dennis", res, jac,
                                                            res_init, jac_init, x_init;
                                                            scale=scale, verbose=verbose,
                                                            print_steps=print_steps)
end
