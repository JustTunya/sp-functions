function f = ContinuousPDF_2D(x, y, distribution_type, ~)
    if ~isequal(size(x), size(y)), error("X must be equal to Y"); end

    f = zeros(size(x));

    switch (distribution_type)
        case "L6"
            mask = (0 <= x & x <= 1) & (-2 <= y & y <= 1);
            f(mask) = (1/11) * ( x(mask) .* (0.5 + 3*y(mask).^2) - (x(mask).^2)/4 + 2 );
        case "L7" 
            mask = (0 <= x & x <= 2) & (-1 <= y & y <= 1);
            f(mask) = 1/4 * ( x(mask) .* y(mask).^2 .* (x(mask).^2 / 4 + 2) + y(mask).^2 / 2);
        otherwise
            error("Unknown type");
    end
end