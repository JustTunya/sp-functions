function F = DiscreteCDF(x, distribution_type, parameters)
    f = @(t) DiscretePDF(t, distribution_type, parameters);

    switch (distribution_type)
        case "bernoulli"
            x_min = 0; x_max = 1;
        case "binomial"
            x_min = 0; x_max = parameters(2);
        case "hypergeometric"
            N = parameters(1); 
            K = parameters(2); 
            n = parameters(3);
            x_min = max(0, n-(N-K));
            x_max = min(n, K);
        case "pascal"
            x_min = 0; x_max = Inf;
        case "poisson"
            x_min = 0; x_max = Inf;
        case "geometric"
            x_min = 1; x_max = Inf;
        case "L2"
            x_min = 1; x_max = 7;
        otherwise
            error("Unknown distribution_type: %s", string(distribution_type));
    end

    wasRow = isrow(x);
    x = x(:);
    if any(abs(x - round(x)) > 0)
        error("x must contain integers only");
    end
    if any(diff(x) < 0)
        error("x must be increasing");
    end

    n = length(x);
    F = zeros(n,1);
 
    if x(1) < x_min
        F(1) = 0;
    else
        hi = min(x(1), x_max);
        if isfinite(hi)
            F(1) = sum(f((x_min:hi).'));
        else
            F(1) = sum(f((x_min:x(1)).'));
        end
    end

    for i = 2:n
        if x(i) < x_min
            F(i) = 0;
        else
            lo = max(x(i - 1) + 1, x_min);
            hi = min(x(i), x_max);
            if lo <= hi
                F(i) = F(i - 1) + sum(f((lo:hi).'));
            else
                F(i) = F(i - 1);
            end
        end
    end

    if wasRow, F = F.'; end
end