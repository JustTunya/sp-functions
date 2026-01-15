function X = CauchyNorm(n)
    X = zeros(1, n);
    accepted = 0;

    while accepted < n
        U = UMersenneTwisterRNG(1);
        V = UMersenneTwisterRNG(1);
        Y = tan(pi * (U - 1/2));
        acc_ratio = (sqrt(exp(1)) / 2) * (1 + Y.^2) * exp(-1/2 * Y.^2);

        if (V <= acc_ratio)
            accepted = accepted + 1;
            X(accepted) = Y;
        end
    end
end 