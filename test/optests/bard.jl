const ybard = [1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1,
               3.5e-1, 3.9e-1, 3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1,
               1.34e0, 2.10e0, 4.39e0]

const bard = let res_init=zeros(15), jac_init=zeros(15,3), x_init=[1.0, 1.0, 1.0]

    function res(x, r)
        for i = 1:15
            u = Float64(i)
            v = 1.6e1 - u
            w = min(u, v)
            r[i] = ybard[i] - (x[1] + u / (x[2] * v + x[3] * w))
        end
        return r
    end

    function jac(x, j)
        for i = 1:15
            j[i, 1] = -1.e0
            u = Float64(i)
            v = 1.6e1 - u
            w = min(u, v)
            t = u/(x[2]*v + x[3]*w)^2
            j[i, 2] = v*t
            j[i, 3] = w*t
        end
        return j
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("bard", res, jac,
                                                   res_init, jac_init, x_init;
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end
