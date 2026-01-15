function f = ContinuousPDF(x, distribution_type, parameters)
    n = length(x);

    switch (distribution_type)
        case "normal"
            mu = parameters(1);
            sigma = parameters(2);

            if (sigma <= 0)
                error("The standard deviation must be a strictly positive number!");
            end

            f = zeros(1,n);
            for i = 1:n
                f(i) = (1 / sqrt(2 * pi) / sigma) * exp(-(x(i) - mu)^2 / (2 * sigma^2));
            end
        case "exponential"
            lambda = parameters(1);
            if lambda <= 0
                error("lambda must be strictly positive!");
            end
            
            f = zeros(1,n);
            for i = 1:n
                if x(i) >= 0
                    f(i) = lambda * exp(-(lambda * x(i)));
                end
            end
        case "beta"
            a = parameters(1);
            b = parameters(2);
            if a <= 0 || b <= 0
                error("a and b must be strictly positive!");
            end
            
            f = zeros(1,n);
            for i = 1:n
                if x(i) > 0 && x(i) < 1
                    f(i) = (x(i)^(a - 1)) * ((1 - x(i))^(b - 1)) / (gamma(a) * gamma(b) / gamma(a + b));
                else
                    f(i) = 0;
                end
            end
        case "student"
            nu = parameters(1);
            if nu <= 0
                error("Degrees of freedom n must be strictly positive!");
            end

            c = gamma((nu + 1) / 2) / (sqrt(nu * pi) * gamma(nu / 2));
            f = zeros(1,n);
            for i = 1:n
                f(i) = c * (1 + (x(i)^2) / nu)^(-(nu + 1) / 2);
            end
        case "chi2"
            r = parameters(1);
            sigma = parameters(2);
            if r < 1 || floor(r) ~= r
                error("r must be a positive integer!");
            end
            if sigma <= 0
                error("sigma must be strictly positive!");
            end
        
            f = zeros(1, n);
            for i = 1:n
                if x(i) > 0
                    num = (x(i) / sigma^2)^(r/2 - 1) * exp(-x(i) / (2*sigma^2));
                    den = 2^(r/2) * gamma(r/2);
                    f(i) = num / den / sigma^2;
                else
                    f(i) = 0;
                end
            end
        case "snedecor-fisher"
            m = parameters(1);
            n2 = parameters(2);
            if m <= 0 || n2 <= 0
                error("m and n must be strictly positive!");
            end
            
            f = zeros(1,n);
            c = (1 / beta(m/2, n2/2)) * (m/n2)^(m/2);
            for i = 1:n
                if x(i) >= 0
                    f(i) = c * (x(i)^(m / 2 - 1)) * (1 + (m / n2)*x(i))^(-(m + n2) / 2);
                else
                    f(i) = 0;
                end
            end
        case "gamma"
            a = parameters(1);
            b = parameters(2);

            if (a <= 0 || b <= 0)
                error("Wrong parameters!")
            end

            f = zeros(1,n);

            for  i = 1 : n
                if (x(i) <= 0)
                    f(i) = 0;
                else
                    f(i) = (1 / ((b^a) * gamma(a))) * (x(i)^(a-1)) * exp(-x(i)/b);
                end
            end
        case "L2"
            alpha = 1/5;

            f = zeros(size(x));
            m1 = (-3 < x) & (x <= -1);
            m2 = (-1 < x) & (x <= 0);
            m3 = (0 < x)  & (x <= 2);
        
            f(m1) = (x(m1).^2 + 1) / 12;
            f(m2) = (x(m2).^2) / 9 .* (1 + x(m2).^2 / 3);
            f(m3) = (alpha / 6) * x(m3);
        case "L6_x"
            alpha = 1/11;

            f = zeros(size(x));
            m = (x >= 0 & x <= 1);
            f(m) = alpha .* (21/2 .* x(m) - 3/4 .* x(m).^2 + 6);
        case "L6_y"
            alpha = 1/11;

            f = zeros(1, n);
            m = (x >= -2 & x <= 1);
            f(m) = alpha .* ( (13/6) + 3/2 .* x(m).^2 );
        case "L12"
            theta = parameters(1);

            if theta <= 0
                error("theta must be strictly positive!");
            end

            f = zeros(1, n);

            for i = 1:n
                if x(i) >= 0 && x(i) <= 2
                    f(i) = (3 * theta / 2) * (1 - x(i) / 2)^(3 * theta - 1);
                else
                    f(i) = 0;
                end
            end
        otherwise
            error("Unknown distribution_type: %s", string(distribution_type));
    end
end