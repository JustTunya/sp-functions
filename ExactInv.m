function X = ExactInv(distribution_type, parameters, n)
    U = UMersenneTwisterRNG(n);
    
    switch (distribution_type)
        case "exponential"
            lambda = parameters(1);
            X = -(1/lambda) * log(1 - U);
        case "cauchy"
            sigma = parameters(1);
            X = sigma * tan(pi * (U - 0.5));
        case "rayleigh"
            sigma = parameters(1);
            X = sigma * sqrt(-2 * log(1 - U));
        case "triangular"
            a = parameters(1);
            X = a * (1 - sqrt(1 - U));
        case "rayleigh_tail"
            a = parameters(1);
            X = sqrt(a^2 - 2 * log(1 - U));
        case "pareto"
            a = parameters(1);
            b = parameters(2);
            X = b * (1 - U).^(-1 / a);
        case "L3"
            u = max(min(U,1),0);
            X = -2 + 0.5 * sqrt(4 + 96*u);
        case "L8"
            theta = parameters(1);
            X = 2 * theta ./ ((1 - U).^(1/5));
        otherwise
            error("Unknown distribution type.");
    end
end