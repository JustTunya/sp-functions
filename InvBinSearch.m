function sample = InvBinSearch(X, uniform_rng, count)
    [row_count, column_count] = size(X);

    if (row_count ~= 2)
        error("Wrong input discrete random variable!");
    end
    if any(X(2,:) < 0 | X(2,:) > 1)
        error("Wrong input discrete random variable!");
    end

    eps = 1e-3;
    total_probability = sum(X(2,:));
    if (total_probability > 1 + eps)
        error("Wrong input discrete random variable!");
    end

    if (total_probability < 1)
        X(2, column_count) = 1 - sum(X(2, 1:column_count - 1));
    end

    q = cumsum(X(2,:));
    n = column_count;

    sample = zeros(1, count);

    for k = 1:count
        switch (uniform_rng)
            case {"LEcuyer", 1}
                u = ULEcuyerRNG();
            case {"MersenneTwister", 2}
                u = UMersenneTwisterRNG();
            otherwise
                u = rand();
        end

        lo = 1; hi = n;
        while lo < hi
            mid = floor((lo + hi) / 2);
            if (u <= q(mid))
                hi = mid;
            else
                lo = mid + 1;
            end
        end
        sample(k) = X(1, lo);
    end
end