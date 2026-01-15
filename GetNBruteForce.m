function n = GetNBruteForce(distribution_type, param, alpha, eps, experiment_count)
    switch distribution_type
        case "normal"
            mu = param(1); sigma = param(2);
            f = @(n) (sigma * BoxMuller(n) + mu);
        case "bernoulli"
            p = param(1); mu = p;
            x = [0, 1; 1 - p, p];
            f = @(n) (InvSeqSearch(x, 'ULEcuyer', n));
        otherwise
            error("Unknown distribution_type: %s", string(distribution_type));
    end

    P = 0; n = 0;

    while (P < 1 - alpha)
        n = n + 1;
        favorable = 0;

        for i = 1 : experiment_count
            X = f(n);

            if (abs(mean(X) - mu) < eps)
                favorable = favorable + 1;
            end
        end

        P = favorable / experiment_count;

        stem(n, P, "m");
        hold on;
        drawnow;
    end
end