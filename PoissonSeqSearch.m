function X = PoissonSeqSearch(params, n)
    lambda = params(1);
    
    X = zeros(1, n);
    U = UMersenneTwisterRNG(n);

    for k = 1 : n
        i = 0;
        p = exp(-lambda); 
        S = p;
        
        while (U(k) > S)
            i = i + 1;
            p = p * (lambda / i);
            S = S + p;
        end
        
        X(k) = i;
    end
end