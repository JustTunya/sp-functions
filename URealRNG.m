function [Y, new_initial_value] = URealRNG(initial_value, generator_type, alpha, beta, n)
    if alpha >= beta
        error("URealRNG:InvalidRange", "alpha must be less than beta");
    end

    switch (generator_type)
        case {"URNG1"}
            m = 2^31 - 1;
            [Y, new_initial_value] = URNG1(initial_value, n);
            Y = Y ./ (m - 1);
        case {"URNG2"}
            m = 2^31 - 1;
            [Y, new_initial_value] = URNG2(initial_value, n);
            Y = Y ./ (m - 1);
        case {"ULEcuyer"}
            n_cell = num2cell(n);
            Y = ULEcuyerRNG(n_cell{1:end}); 
            new_initial_value = [];
        case {"UMersenneTwister"}
            n_cell = num2cell(n);
            Y = UMersenneTwisterRNG(n_cell{1:end});
            new_initial_value = [];
        otherwise
            error("Unknown generator type: %s", generator_type);
    end
    
    Y = alpha + Y .* (beta - alpha);
end
