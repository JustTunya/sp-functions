function X = LaplaceNorm(n, mu, sigma)
    if nargin < 2 || isempty(mu)
        mu = 1;
    end
    if nargin < 3 || isempty(sigma)
        sigma = 1;
    end

    X = zeros(1, n);

    for i = 1 : n
        Y = 0; V = 1;
        
        while (Y - 1)^2 > -2 * log(abs(V))
            Y = ExactInv("exponential", 1, 1);
            V = UMersenneTwisterRNG(1) * 2 - 1;
        end

        X(i) = mu + sigma * Y * sign(V);
    end
end