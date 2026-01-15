function [Y1, Y2] = CorrelatedNorm2D(mu1, mu2, sigma1, sigma2, p, n)
    if (sigma1 <= 0)
        error("sigma1 should be positive");
    end
    if (sigma2 <= 0)
        error("sigma2 should be positive");
    end
    if (abs(p) > 1)
        error("p must be in [-1, 1]");
    end

    mu = [mu1; mu2];
    L = [sigma1, 0; p * sigma2, sigma2 * sqrt(max(0, 1 - p^2))];
    Y = zeros(2, n);

    for k = 1:n
        U1 = 1 - ULEcuyerRNG(1);
        U2 = ULEcuyerRNG(1);

        R = sqrt(-2 * log(U1));
        Theta = 2 * pi * U2;
        X = [R * cos(Theta); R * sin(Theta)];
        Y(:,k) = mu + L * X;
    end

    Y1 = Y(1,:);
    Y2 = Y(2,:);
end