function [ci_u, ci_mu, u_value, p_value, H] = UTest(X, mu_0, sigma, alpha, tail)
    if ~(sigma > 0)
        error('Sigmas should be positive!');
    end
    if ~(0 < alpha && alpha < 1)
        error('Alpha should be in the interval (0, 1).');
    end

    n = length(X); 
    s = sigma / sqrt(n);
    X_ = mean(X);
    u_value = (X_ - mu_0) / s;

    switch (tail)
        case {'both', 0}
            u = norminv(1 - alpha / 2, 0, 1);
            ci_u(1) = -u; ci_u(2) =  u;
            ci_mu(1) = X_ - s * u; ci_mu(2) = X_ + s * u;
            
            p_value = 2.0 * normcdf(-abs(u_value), 0, 1);
        case {'right', 1}
            u = norminv(1 - alpha, 0, 1);
            ci_u(1) = -inf; ci_u(2) =  u;
            ci_mu(1) = X_ - s * u; ci_mu(2) = inf;
            
            p_value = 1.0 - normcdf(u_value, 0, 1);
        case {'left', -1}
            u = norminv(alpha, 0, 1);
            ci_u(1) = u; ci_u(2) = inf;
            ci_mu(1) = -inf; ci_mu(2) = X_ - s * u;
            
            p_value = normcdf(u_value, 0, 1);
        otherwise
            error("Unknown tail!");
    end

    H = ~(u_value > ci_u(1) && u_value < ci_u(2));

    if (p_value < alpha)
        disp('Warning: small p-value cast doubt on the validity of the null-hypothesis!');
    end
end