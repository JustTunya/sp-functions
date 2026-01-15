function [X, new_initial_value] = URNG2(initial_value, n)
    m = 2^31 - 1; a = 2^16 + 3; c = 0;
    initial_value = mod(initial_value, m);
    if initial_value <= 0
        initial_value = 1;
    end
    if mod(initial_value, 2) == 0
        initial_value = initial_value + 1;
    end
    [X, new_initial_value] = LinearCongruentialGenerator(m, a, c, initial_value, n);
    if mod(new_initial_value, 2) == 0
        new_initial_value = new_initial_value + 1;
    end
end