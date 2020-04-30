const ymeyer = [3.478e4, 2.861e4, 2.365e4, 1.963e4, 1.637e4, 1.372e4,
                1.154e4, 9.744e3, 8.261e3, 7.030e3, 6.005e3, 5.147e3,
                4.427e3, 3.820e3, 3.307e3, 2.872e3]

const meyer = let res_init=zeros(16), jac_init=zeros(16,3), x_init=[0.02, 4.e3, 2.5e2]

    function res(x, r)
        for i = 1:16
            ti = Float64(5*i + 45)
            r[i] = x[1] * exp(x[2] / (ti + x[3])) - ymeyer[i]
        end
        return r
    end

    function jac(x, j)
        for i = 1:16
            ti = Float64(5*i + 45)
            u = ti + x[3]
            t = exp(x[2] / u)
            j[i, 1] = t
            j[i, 2] = x[1] * t/u
            j[i, 3] = -x[1] * x[2] * t / (u * u)
        end
        return j
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("meyer", res, jac,
                                                   res_init, jac_init, x_init;
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end
