function f = DiscretePDF(x, distribution_type, parameters)
    sort(x);
    x = round(x);
    n = length(x);

    switch (distribution_type)
        case "bernoulli"
            p = parameters(1);
            if p < 0 || p > 1
                error("p must be between 0 and 1");
            end
            
            f = zeros(1,n);
            for i = 1:n
                if x(i) == 1
                    f(i) = p;
                elseif x(i) == 0
                    f(i) = 1 - p;
                else
                    f(i) = 0;
                end
            end
        case "binomial"
            p = parameters(1);
            r = parameters(2);

            f = zeros(1, n);
            for i=1:n
                if(x(i) < 0)
                    error("Incorrect input data");
                end
                k = x(i);
                f(i) = nchoosek(r, k) * (p^k) * (1 - p)^(r - k);
            end
        case "hypergeometric"
            N = parameters(1); K = parameters(2); ntr = parameters(3);
            if (N < 1 || K < 0 || K > N || ntr < 0 || ntr > N)
                error("Wrong parameters");
            end
            
            f = zeros(1, n);
            for i = 1:n
                if (max(0, ntr - N + K) <= x(i) && x(i) <= min(ntr, K))
                    f(i) = (nchoosek(K, x(i)) * nchoosek(N - K, ntr - x(i))) / nchoosek(N, ntr);
                end
            end
        case "pascal"
            r = parameters(1);
            p = parameters(2);
            if (r <= 1 || p < 0 || p > 1)
                error("Wrong parameters");
            end
            
            f = zeros(1, n);
            for i = 1:n
                if x(i) < 0
                    error("Wrong inputs");
                end
                f(i) = nchoosek(r + x(i) - 1, x(i)) * (p^r) * (1 - p)^(x(i));
            end
        case "poisson"
            lambda = parameters(1);
            if lambda <= 0
                error("lambda must be positive");
            end
            
            f = zeros(1,n);
            for i = 1:n
                k = x(i);
                if k >= 0
                    f(i) = exp(-lambda) * (lambda^k) / factorial(k);
                else
                    f(i) = 0;
                end
            end
        case "geometric"
            p = parameters(1);
            
            if (p < 0 || p > 1)
                error("Wrong parameter");
            end
            
            f = zeros(1, n);

            q = 1 - p;
            for i = 1:n
                if (x(i) < 1)
                    error("Incorrect input data");
                else
                    f(i) = q^(x(i) - 1) * p;
                end
            end
        case "L2"
            x_vals = 1:7;
            raw = abs(sin(x_vals * pi / 8)) + (x_vals / 30);
            probs = raw / sum(raw);

            f = zeros(size(x));
            for i = 1:length(x)
                idx = find(x_vals == x(i), 1);
                if ~isempty(idx)
                    f(i) = probs(idx);
                else
                    f(i) = 0;
                end
            end
        otherwise
            error("Unknown distribution_type: %s", string(distribution_type));
    end
end