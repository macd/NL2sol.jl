const davidon = let res_init=zeros(15), jac_init=zeros(15,15), x_init=fill(0.1, 15)

    function res(x, r)
        for i = 1:14
            r1 = 0.0
            ti = Float64(i)
            t = 1.0
            for j = 1:15
                r1 += t * x[j]
                t *= ti
            end
            r[i] = r1
        end
        r[15] = x[1] - 1.0
        return r
    end

    function jac(x, jac)
        for i = 1:14
            ti = Float64(i)
            t = 1.0
            for k = 1:15
                jac[i, k] = t
                t *= ti
            end
        end
        jac[15, 1] = 1.0
        for k = 2:15
            jac[15, k] = 0.0
        end
        return jac
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("davidon", res, jac,
                                                   res_init, jac_init, x_init;
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end
