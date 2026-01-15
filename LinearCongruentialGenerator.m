function [X, new_initial_value] = LinearCongruentialGenerator(m, a, c, initial_value, n)
    X = zeros(1, n);
    X(1) = initial_value;

    for i = 2:n
        X(i) = mod(a * X(i-1) + c, m);
    end

    new_initial_value = mod(a * X(n) + c, m);
end