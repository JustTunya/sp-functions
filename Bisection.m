function X = Bisection(distribution_type, parameters, delta, n, x_min, x_max)
    F = @(x) ContinuousCDF(x, distribution_type, parameters);
    U = UMersenneTwisterRNG(n);
    X = zeros(1, n);
    
    for i = 1 : n
        a = x_min;
        b = x_max;
        x = (a + b) / 2;
        while (b - a > delta && abs(U(i) - F(x)) > delta)
            if (U(i) < F(x))
                b = x;
            else
                a = x;
            end
            x = (a + b) / 2;
        end
        X(i) = x;
    end
end