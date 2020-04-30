const jennrich_sampson = let res_init=zeros(10), jac_init=zeros(10,2), x_init=[0.3, 0.4]

    function res(x, r)
        for i = 1:10
            ti = Float64(i)
            r[i] = 2.0e0 + 2.0e0*ti - (exp(ti*x[1]) + exp(ti*x[2]))
        end
        return r
    end

    function jac(x, j)
        for i = 1:10
            ti = Float64(i)
            j[i,1] = -ti * exp(ti*x[1])
            j[i,2] = -ti * exp(ti*x[2])
        end
        return j
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("jennrich_sampson", res, jac,
                                                            res_init, jac_init, x_init;
                                                            scale=scale, verbose=verbose,
                                                            print_steps=print_steps)
end
