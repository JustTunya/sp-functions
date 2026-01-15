function F = ContinuousCDF(x, distribution_type, parameters)
    f = @(t) ContinuousPDF(t, distribution_type, parameters);

    switch distribution_type
        case {"normal","student"}
            x_min = -Inf;
        case {"exponential","beta","chi2","snedecor-fisher","gamma","L6_x","L12"}
            x_min = 0;
        case "L6_y"
            x_min = -2;
        case "L2"
            x_min = -3;
        otherwise
            error("Unknown distribution_type: %s", string(distribution_type));
    end

    if ismember(distribution_type, ["L6_x","L6_y"])
        x = max(x, x_min);
    end

    n = numel(x); F = zeros(1,n);

    if (n == 0)
        return;
    end

    F(1) = integral(f, x_min, x(1));

    for i=2:n
        F(i) = F(i-1) + integral(f, x(i-1), x(i));
    end
end