function [ci_chi2, ci_std, chi2_value, p_value, H] = Chi2Test(X, sigma_0, alpha, tail)
    if (~isnumeric(X))
        error('The input samples X must be a numeric array!');
    end
    if (~isvector(X))
        error('The input samples X must be a one-dimensional array (vector)!');
    end

    n = length(X);
    sigma_2_ = var(X);
    chi2_value = (n-1) * sigma_2_ / sigma_0^2;

    switch (tail)
        case {'both', 0}
            ci_chi2(1) =  chi2inv(alpha/2, n-1); ci_chi2(2) =  chi2inv(1 - alpha/2, n-1);
            ci_std(1) = sqrt((n - 1) * sigma_2_ / ci_chi2(2)); ci_std(2) = sqrt((n - 1) * sigma_2_ / ci_chi2(1));
            
            probability = chi2cdf(chi2_value, n-1);
            p_value = 2.0 * min(probability, 1 - probability);
        case {'right', 1}
            ci_chi2(1) = 0; ci_chi2(2) = chi2inv(1 - alpha, n-1);   
            ci_std(1) = sqrt((n - 1) * sigma_2_ / ci_chi2(2)); ci_std(2) = inf;
                
            probability = chi2cdf(chi2_value, n-1);
            p_value = 1 - probability; 
        case {'left', -1}
            ci_chi2(1) = chi2inv(alpha, n-1); ci_chi2(2) = inf;
            ci_std(1) = 0; ci_std(2) = sqrt((n - 1) * sigma_2_ / ci_chi2(1));
                
            probability = chi2cdf(chi2_value, n-1);
            p_value = probability;
        otherwise
            error("Unknown tail!");
    end

    H = ~(chi2_value > ci_chi2(1) && chi2_value < ci_chi2(2));

    if (p_value < alpha)
        disp('Warning: small p-value cast doubt on the validity of the null-hypothesis!');
    end
end