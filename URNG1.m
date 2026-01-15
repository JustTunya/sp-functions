function [X, new_initial_value] = URNG1(initial_value, n)
  m = 2^31-1; a = 7^5; c = 0;
  [X, new_initial_value] = LinearCongruentialGenerator(m, a, c, initial_value, n);
end