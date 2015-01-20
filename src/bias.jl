#collection of bias correction routines

ma(K,N) = ma(K,N,2.0)
function ma(K, N, base)
    return (K-1)/(2*N*log(base))
end

blank(args...) = 0
