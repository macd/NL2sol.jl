using NL2sol

function testone(prb, nlres, nljac, r_init, j_init, x_init; scale=1,
                 mxfc=200, mxiter=200, verbose=false, print_steps=false)
    n = length(r_init)
    p = length(x_init)
    iv, v = nl2_set_defaults(n, p)
    iv[MXFCAL] = mxfc
    iv[MXITER] = mxiter
    iv[PRUNIT] = (verbose && print_steps) ? 6 : 0
    xscale = [1.0, 10.0, 100.0]
    pass = true
    for i = 1:scale
        verbose && println("Running NL2sol on problem ", prb, " at scale ", xscale[i])
        results = nl2sol(nlres, nljac, xscale[i] * x_init, n, iv, v)
        verbose && println(return_code[iv[1]])
        verbose && println(results)
        pass = pass && iv[1] in (3, 4, 5, 6)
    end
    return pass
end
