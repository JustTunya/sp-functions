function X = LaplaceNorm(n)
    X = zeros(1, n);

    for i = 1 : n
        Y = 0; V = 1;
        
        while (Y - 1)^2 > -2 * log(abs(V))
            Y = ExactInversion("exponential", 1, 1);
            V = UMersenneTwisterRNG(1) * 2 - 1;
        end

        X(i) = Y * sign(V);
    end
end