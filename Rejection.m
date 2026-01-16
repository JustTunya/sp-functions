function [X, acceptance_rate] = Rejection(n, f, limits, M)
    X = zeros(1, n);
    i = 1; tries = 0;

    while i <= n
        uv = ULEcuyerRNG(1, 2);
        U = uv(1); V = uv(2);
        Y = limits(1) + (limits(2) - limits(1)) * V;
        tries = tries + 1;

        if (U * M <= f(Y))
            X(i) = Y;
            i = i + 1;
        end
    end

    acceptance_rate = n / tries;
end
