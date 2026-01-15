function [ci_t, ci_mu, t_value, p_value, H] = TTest(X, mu_0, alpha, tail)
    if ~isnumeric(X)
        error('The sample X must be numeric!');
    end
    if ~isvector(X)
        error('The sample X must be a one-dimensional array (vector)!');
    end

    X = X(:);
    n = length(X);
    X_ = sum(X) / n;
    sigma_ = sqrt(sum((X - X_).^2) / (n-1));
    s = sigma_ / sqrt(n);
    t_value = (X_ - mu_0) / s;

    switch (tail)
        case {'both', 0}
            t = tinv(1 - alpha/2, n-1);
            ci_t(1) = -t; ci_t(2) =  t;
            ci_mu(1) = X_ - s * t; ci_mu(2) = X_ + s * t;
            
            p_value = 2.0 * tcdf(-abs(t_value), n-1);
        case {'right', 1}
            t = tinv(1 - alpha, n - 1);
            ci_t(1) = -inf; ci_t(2) =  t;
            ci_mu(1) = X_ - s * t; ci_mu(2) = inf;
            
            p_value = 1.0 - tcdf(t_value, n - 1);                 
        case {'left', -1}
            t = tinv(alpha, n - 1);
            ci_t(1) = t; ci_t(2) = inf;
            ci_mu(1) = -inf; ci_mu(2) = X_ - s * t;
            
            p_value = tcdf(t_value, n - 1);
        otherwise
            error("Unknown tail!");
    end

    H = ~(t_value > ci_t(1) && t_value < ci_t(2));

    if (p_value < alpha)
        disp('Warning: small p-value cast doubt on the validity of the null-hypothesis!');
    end
end