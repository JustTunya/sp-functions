function F = ContinuousCDF_2D(x, y, distribution_type, params)
    if ~isequal(size(x), size(y)), error("X must be equal to Y"); end

    f = @(u, v) ContinuousPDF_2D(u, v, distribution_type, params);
    F = zeros(size(x));

    switch (distribution_type)
        case "L6"
            u_min = 0; u_max = 1;
            v_min = -2; v_max = 1;

            for i = 1:numel(x)
                xi = x(i); yi = y(i);

                if (xi <= u_min || yi <= v_min)
                    F(i) = 0;
                else
                    u_upper = min(xi, u_max);
                    v_upper = min(yi, v_max);

                    F(i) = integral2(f, u_min, u_upper, v_min, v_upper);
                end
            end
        otherwise
            error("Unknown distribution_type: %s", string(distribution_type));
    end
end