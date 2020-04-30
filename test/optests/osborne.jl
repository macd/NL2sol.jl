const yosb1 = [8.44e-1, 9.08e-1, 9.32e-1, 9.36e-1, 9.25e-1, 9.08e-1,
               8.81e-1, 8.50e-1, 8.18e-1, 7.84e-1, 7.51e-1, 7.18e-1,
               6.85e-1, 6.58e-1, 6.28e-1, 6.03e-1, 5.80e-1, 5.58e-1,
               5.38e-1, 5.22e-1, 5.06e-1, 4.90e-1, 4.78e-1, 4.67e-1,
               4.57e-1, 4.48e-1, 4.38e-1, 4.31e-1, 4.24e-1, 4.20e-1,
               4.14e-1, 4.11e-1, 4.06e-1]

const osborne1 = let res_init=zeros(33), jac_init=zeros(33,5), x_init=[.5, 1.5, -1.0, .01, .02]

    function res(x, r)
        for i = 1:33
            ti = 1.0e1 * Float64(1-i)
            r[i] = yosb1[i] - (x[1] + x[2] * exp(x[4] * ti) + x[3] * exp(x[5] * ti))
        end
        return r
    end

    function jac(x, j)
        for i = 1:33
            ti = 1.0e1 * Float64(1-i)
            j[i, 1] = -1.e0
            j[i, 2] = -exp(x[4] * ti)
            j[i, 3] = -exp(x[5] * ti)
            j[i, 4] = ti * x[2] * j[i, 2]
            j[i, 5] = ti * x[3] * j[i, 3]
        end
        return j
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("osborne", res, jac,
                                                            res_init, jac_init, x_init;
                                                            scale=scale, verbose=verbose,
                                                            print_steps=print_steps)
end


const yosb2 = [ 1.366e0, 1.191e0, 1.112e0, 1.013e0, 9.91e-1, 8.85e-1,
                8.31e-1, 8.47e-1, 7.86e-1, 7.25e-1, 7.46e-1, 6.79e-1,
                6.08e-1, 6.55e-1, 6.16e-1, 6.06e-1, 6.02e-1, 6.26e-1,
                6.51e-1, 7.24e-1, 6.49e-1, 6.49e-1, 6.94e-1, 6.44e-1,
                6.24e-1, 6.61e-1, 6.12e-1, 5.58e-1, 5.33e-1, 4.95e-1,
                5.00e-1, 4.23e-1, 3.95e-1, 3.75e-1, 3.72e-1, 3.91e-1,
                3.96e-1, 4.05e-1, 4.28e-1, 4.29e-1, 5.23e-1, 5.62e-1,
                6.07e-1, 6.53e-1, 6.72e-1, 7.08e-1, 6.33e-1, 6.68e-1,
                6.45e-1, 6.32e-1, 5.91e-1, 5.59e-1, 5.97e-1, 6.25e-1,
                7.39e-1, 7.10e-1, 7.29e-1, 7.20e-1, 6.36e-1, 5.81e-1,
                4.28e-1, 2.92e-1, 1.62e-1, 9.8e-2, 5.4e-2]

# uftolg is a machine-dependent constant.  it is just slightly
# larger than the log of the smallest positive machine number.
const uftolg = 1.999e0 * log(sqrt(1.001 * 2.2250738585072014e-308))

const osborne2 = let res_init=zeros(65), jac_init=zeros(65,11),
                     x_init=[1.3, .65, .65, .7, .6, 3., 5., 7., 2., 4.5, 5.5]

    function res(x, r)
        for i = 1:65
            ti = 0.1e0 * Float64(1-i)
            ri = x[1] * exp(x[5] * ti)
            for j = 2:4
                t = 0.e0
                theta = -x[j+4] * (ti + x[j+7]) ^ 2
                if theta > uftolg
                    t = exp(theta)
                end
                ri = ri + x[j] * t
            end
            r[i] = yosb2[i] - ri
        end
        return r
    end

    function jac(x, j)
        for i = 1:65
            ti = Float64(1 - i) * 1.e-1
            j[i, 1] = -exp(x[5] * ti)
            j[i, 5] = x[1] * ti * j[i,1]
            for k = 2:4
                t = x[k + 7] + ti
                theta = -x[k+4]*t*t
                r2 = theta > uftolg ? -exp(theta) : 0.e0
                j[i, k] = r2
                r2 = -t * r2 * x[k]
                j[i, k+4] = r2*t
                j[i, k+7] = 2.e0 * x[k+4] * r2
            end
        end
        return j
    end

    f(;scale=1, verbose=false, print_steps=false) = testone("osborne2", res, jac,
                                                            res_init, jac_init, x_init;
                                                            scale=scale, verbose=verbose,
                                                            print_steps=print_steps)
end
