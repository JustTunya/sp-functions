function [p_min, p_max] = GetSafetyInterval(distribution_type, X, alpha)
    x = norminv(1 - alpha/2, 0, 1);
    n = length(X);
    X_mean = mean(X);
    
    switch distribution_type
        case "binomial"
            a = 64 * n - 8 * x^2;
            b = -16 * n * X_mean - 8 * x^2;
            c = n * X_mean^2;
        case "poisson"
            a = n;
            b = -(2 * n * X_mean + x^2);
            c = n * X_mean^2;
        case "L8"
            a = (25/4) - ((5 * x^2) / (12 * n));
            b = -5 * X_mean; 
            c = X_mean^2;
        otherwise
            error("Unknown distribution type!");
    end
    
    if (a <= 0)
        error("a is not positive.");
    end
    
    delta = b^2 - 4*a*c;
    
    if (delta < 0)
        error("delta is negative. No real roots.");
    end
    
    roots = [(-b + sqrt(delta)) / (2*a), (-b - sqrt(delta)) / (2*a)];
    roots = sort(roots);
    
    p_min = roots(1); p_max = roots(2);
end