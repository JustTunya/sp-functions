function [ci_f, ci_lambda, f_value, p_value, H] = FTest2D(X, Y, alpha, tail)
    assert(0 < alpha && alpha < 1, "alpha must be between 0 and 1.");
    
    m = length(X); n = length(Y);
    f_value = var(X) / var(Y);
    lambda_X_Y_ = sqrt(f_value);
    
    switch (tail)
        case {'both', 0}
            ci_f(1) = finv(alpha/2, m-1, n-1); ci_f(2) = finv(1 - alpha/2, m-1, n-1);
            ci_lambda(1) = (1 / sqrt(ci_f(2))) * lambda_X_Y_; ci_lambda(2) = (1 / sqrt(ci_f(1))) * lambda_X_Y_;
              
            probability = fcdf(f_value, m-1, n-1);
            p_value = 2.0 * min(probability, 1 - probability);
        case {'right', 1}
            ci_f(1) = 0; ci_f(2) = finv(1 - alpha, m-1, n-1);
            ci_lambda(1) = (1 / sqrt(ci_f(2))) * lambda_X_Y_; ci_lambda(2) = inf;
              
            probability = fcdf(f_value, m-1, n-1);
            p_value = 1 - probability;
        case {'left', -1}
            ci_f(1) = finv(alpha, m-1, n-1); ci_f(2) = inf;
            ci_lambda(1) = 0; ci_lambda(2) = (1 / sqrt(ci_f(1))) * lambda_X_Y_;
              
            probability = fcdf(f_value, m-1, n-1);
            p_value = probability;
        otherwise
            error("Unknown tail!");
    end
    
    H = ~(f_value > ci_f(1) && f_value < ci_f(2));
    
    if (p_value < alpha)
        warning('Small p-value cast doubt on the validity of the null-hypothesis!');
    end
end