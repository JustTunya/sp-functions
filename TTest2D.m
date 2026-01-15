function [ci_t, ci_delta, t_value, p_value, H] = TTest2D(X, Y, equal_std, alpha, tail)
    if (~isnumeric(X) || ~isnumeric(Y))
        error('The input samples X and Y must be numeric arrays!');
    end
    if (~isvector(X) || ~isvector(Y))
        error('The input samples X and Y must be one-dimensional arrays (vectors)!');
    end

    X = X(:); Y = Y(:);
    m = length(X); n = length(Y);
    sigma_1_2_ = var(X); sigma_2_2_ = var(Y);

    if (equal_std)
        s = sqrt(((m-1) * sigma_1_2_ + (n-1) * sigma_2_2_) * (m + n) / m / n / (m + n - 2));
    else
        s = sqrt(sigma_1_2_ / m + sigma_2_2_ / n);
    end

    if (equal_std)
        eta = m + n - 2;
    else
        w = n * sigma_1_2_ / (n * sigma_1_2_ + m * sigma_2_2_);
        eta = (m - 1) * (n - 1) / ((m - 1) * (1 - w)^2 + (n - 1) * w^2);
    end

    X_Y_ = mean(X) - mean(Y);
    t_value = X_Y_ / s;

    switch (tail)
        case {'both', 0}
            t = tinv(1 - alpha / 2, eta);
            ci_t(1) = -t; ci_t(2) =  t;
            ci_delta(1) = X_Y_ - s * t; ci_delta(2) = X_Y_ + s * t;
            
            p_value = 2.0 * tcdf(-abs(t_value), eta);
        case {'right', 1}
            t = tinv(1 - alpha, eta);
            ci_t(1) = -inf; ci_t(2) =  t;
            ci_delta(1) = X_Y_ - s * t; ci_delta(2) = inf;
        
            p_value = 1.0 - tcdf(t_value, eta);
        case {'left', -1}
            t = tinv(alpha, eta);
            ci_t(1) = t; ci_t(2) = inf;
            ci_delta(1) = -inf; ci_delta(2) = X_Y_ - s * t;
        
            p_value = tcdf(t_value, eta);
        otherwise
            error("Unknown tail!");
    end

    H = ~(t_value > ci_t(1) && t_value < ci_t(2));

    if (p_value < alpha)
        disp('Warning: small p-value cast doubt on the validity of the null-hypothesis!');
    end
end