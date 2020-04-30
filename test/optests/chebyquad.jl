const chebyquad = let res_init=zeros(8), jac_init=zeros(8,8), x_init=[1:8;]/9

    function res(x, r)
        n = length(r)
        r[:] .= 0.0
        for j = 1:n
            tim1 = 1.0e0
            ti = 2.0e0 * x[j] - 1.0e0
            z = ti + ti
            for i = 1:n
                r[i] = r[i] + ti
                tip1 = z * ti - tim1
                tim1 = ti
                ti = tip1
            end
        end
        floatn = Float64(n)
        for i = 1:n
            ti = 0.0e0
            if mod(i,2) == 0
                ti = -1.0e0 / Float64(i*i - 1)
            end
            r[i] = ti - r[i]/floatn
        end
        return r
    end

    function jac(x, j)
        n = size(j)[1]
        for k = 1:n
            tim1 = -1.e0 / Float64(n)
            z = 2.e0 * x[k] - 1.e0
            ti = z * tim1
            tpim1 = 0.e0
            tpi = 2.e0 * tim1
            z = z + z
            for i = 1:n
                j[i, k] = tpi
                tpip1 = 4.e0 * ti + z * tpi - tpim1
                tpim1 = tpi
                tpi = tpip1
                tip1 = z * ti - tim1
                tim1 = ti
                ti = tip1
            end
        end
        return j
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("chebyquad", res, jac,
                                                   res_init, jac_init, x_init;
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end
