function [ci_u, ci_delta, u_value, p_value, H] = UTest2D(X, Y, sigma_1, sigma_2, alpha, tail)
    if ~(sigma_1 > 0 && sigma_2 > 0)
        error('Sigmas should be positive!');
    end
    if ~(0 < alpha && alpha < 1)
        error('Alpha should be in the interval (0, 1).');
    end

    m = length(X); n = length(Y);
    s = sqrt(sigma_1^2 / m + sigma_2^2 / n);
    X_Y_ = mean(X) - mean(Y);
    u_value = X_Y_ / s;

    switch (tail)
        case {'both', 0}
            u = norminv(1 - alpha / 2, 0, 1);
            ci_u(1) = -u; ci_u(2) =  u;
            ci_delta(1) = X_Y_ - s * u; ci_delta(2) = X_Y_ + s * u;
            
            p_value = 2.0 * normcdf(-abs(u_value), 0, 1);
        case {'right', 1}
            u = norminv(1 - alpha, 0, 1);
            ci_u(1) = -inf; ci_u(2) = u;
            ci_delta(1) = X_Y_ - s * u; ci_delta(2) = inf;

            p_value = 1.0 - normcdf(u_value, 0, 1);
        case {'left', -1}
            u = norminv(alpha, 0, 1);
            ci_u(1) = u; ci_u(2) = inf;
            ci_delta(1) = -inf; ci_delta(2) = X_Y_ - s * u;

            p_value = normcdf(u_value,0,1);
        otherwise
            error("Unknown tail!");
    end

    H = ~(u_value > ci_u(1) && u_value < ci_u(2));

    if (p_value < alpha)
        disp('Warning: small p-value cast doubt on the validity of the null-hypothesis!');
    end
end