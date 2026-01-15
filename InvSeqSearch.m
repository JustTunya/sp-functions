function sample = InvSeqSearch(X, uniform_rng, count)
    [row_count, column_count] = size(X);
    
    if (row_count ~= 2)
        error("Wrong input discrete random variable!");
    end
    
    for j = 1:column_count
        if (X(2,j) < 0 || X(2,j) > 1)
            error("Wrong input discrete random variable!");
        end
    end
    
    total_probability = sum(X(2,:));
    if (total_probability > 1.0001)
        error("Wrong input discrete random variable!");
    end
    
    if (total_probability < 1)
        X(2, column_count) = 1 - sum(X(2, 1:column_count-1));
    end
    
    cumulative_probabilities = cumsum(X(2,:));
    sample = zeros(1, count);
    
    for k = 1:count
        switch (uniform_rng)
            case {'LEcuyer', 1}
                u = ULEcuyerRNG();
    
            case {'MersenneTwister', 2}
                u = UMersenneTwisterRNG();
    
            otherwise
                u = rand();
        end
    
        i = 1;
        while (u > cumulative_probabilities(i))
            i = i + 1;
        end
        
        sample(k) = X(1, i);
    end
end