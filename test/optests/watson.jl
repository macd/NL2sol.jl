function watson_res(x, r, p)
    # for watson nex = 13|14|15|16
    for i = 1:29
        ti = Float64(i) / 2.9e1
        r1 = 0.0e0
        r2 = x[1]
        t = 1.0e0
        for j = 2:p
              r1 = r1 + Float64(j-1) * t * x[j]
              t = t * ti
              r2 = r2 + t * x[j]
        end
        r[i] = r1 - r2 * r2 - 1.0e0
    end
    r[30] = x[1]
    r[31] = x[2] - x[1] ^ 2 - 1.0e0
    return r
end


function watson_jac(x, j, p)
    for i = 1:29
        ti = Float64(i) / 2.9e1
        r2 = x[1]
        t = 1.e0
        for k = 2:p
             t = t*ti
             r2 = r2 + t*x[k]
        end
        r2 = -2.e0 * r2
        j[i, 1] = r2
        t = 1.e0
        r2 = ti * r2
        for k = 2:p
            j[i, k] = t * (Float64(k-1) + r2)
            t = t * ti
        end
    end

    for i = 30:31
        for k = 2:p
            j[i, k] = 0.e0
        end
    end
    j[30, 1] = 1.e0
    j[31, 1] = -2.e0*x[1]
    j[31, 2] = 1.e0
    return j
end


const watson6 = let res_init=zeros(31), jac_init=zeros(31,6), x_init=zeros(6)

    res(x, r) = watson_res(x, r, 6)
    jac(x, j) = watson_jac(x, j, 6)

    f(;scale=1, verbose=false, print_steps=false) = testone("watson6", res, jac,
                                                   res_init, jac_init, x_init;
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end

const watson9 = let res_init=zeros(31), jac_init=zeros(31,9), x_init=zeros(9)

    res(x, r) = watson_res(x, r, 9)
    jac(x, j) = watson_jac(x, j, 9)

    f(;scale=1, verbose=false, print_steps=false) = testone("watson9", res, jac,
                                                   res_init, jac_init, x_init;
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end

const watson12 = let res_init=zeros(31), jac_init=zeros(31,12), x_init=zeros(12)

    res(x, r) = watson_res(x, r, 12)
    jac(x, j) = watson_jac(x, j, 12)

    f(;scale=1, verbose=false, print_steps=false) = testone("watson12", res, jac,
                                                   res_init, jac_init, x_init;
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end


# the Fortran driver and test code places a 15 iteration limit on 
# watson20 so we will do so as well to match.
const watson20 = let res_init=zeros(31), jac_init=zeros(31,20), x_init=zeros(20)

    res(x, r) = watson_res(x, r, 20)
    jac(x, j) = watson_jac(x, j, 20)

    f(;scale=1, verbose=false, print_steps=false) = testone("watson20", res, jac,
                                                   res_init, jac_init, x_init;
                                                   #mxiter=15, scale=scale, verbose=verbose,
                                                   scale=scale, verbose=verbose,
                                                   print_steps=print_steps)
end
