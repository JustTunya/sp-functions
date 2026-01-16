function X = GeoTrunc(p, n)
    if (n < 1 || p <= 0 || p >= 1)
        error("Wrong parameters!")
    end

    X = zeros(1,n);
    for i = 1:n
        % U = ULEcuyerRNG();
        % X(i) = log(1-U) / log(1-p);
        % X(i) = -1/(1-exp(-p))*log(U);
    
        lambda = -log(1-p);
        X(i) = ExactInv("exponential", lambda, 1);
        X(i) = ceil(X(i));
    end
end