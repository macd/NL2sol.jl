function brown_res(x, r)
    n = length(x)
    t = x[1] - Float64(n + 1)
    for i = 2:n
        t = t + x[i]
    end
    nm1 = n - 1
    for i = 1:nm1
        r[i] = t + x[i]
    end
    t = x[1]
    for i = 2:n
        t = t * x[i]
    end
    r[n] = t - 1.0e0
    return r
end

# note for these problems m == p
function brown_jac(x, j, n)
    nm1 = n - 1
    for k = 1:n
        for i = 1:nm1
            j[i, k] = 1.0e0
            if i == k
                j[i, k] = 2.0e0
            end
        end
    end
    for k = 1:n
        t = 1.0e0
        for i = 1:n
            t =  i != k ? t * x[i] : t
        end
        j[n, k] = t
    end
    return j
end


const brown5 = let res_init=zeros(5), jac_init=zeros(5,5), x_init=fill(0.5,5)
    jac(x, jac) = brown_jac(x, jac, 5)
    f(;scale=1, verbose=false, print_steps=false) = testone("brown5", brown_res, jac,
                                                            res_init, jac_init, x_init;
                                                            scale=scale, verbose=verbose,
                                                            print_steps=print_steps)
end

const brown10 = let res_init=zeros(10), jac_init=zeros(10,10), x_init=fill(0.5,10)
    jac(x, jac) = brown_jac(x, jac, 10)
    f(;scale=1, verbose=false, print_steps=false) = testone("brown10", brown_res, jac,
                                                            res_init, jac_init, x_init;
                                                            scale=scale, verbose=verbose,
                                                            print_steps=print_steps)
end

const brown30 = let res_init=zeros(30), jac_init=zeros(30,30), x_init=fill(0.5, 30)
    jac(x, jac) = brown_jac(x, jac, 30)
    f(;scale=1, verbose=false, print_steps=false) = testone("brown30", brown_res, jac,
                                                            res_init, jac_init, x_init;
                                                            scale=scale, verbose=verbose,
                                                            print_steps=print_steps)
end

const brown40 = let res_init=zeros(40), jac_init=zeros(40,40), x_init=fill(0.5, 40)
    jac(x, jac) = brown_jac(x, jac, 40)
    f(;scale=1, verbose=false, print_steps=false) = testone("brown40", brown_res, jac,
                                                            res_init, jac_init, x_init;
                                                            scale=scale, verbose=verbose,
                                                            print_steps=print_steps)
end


