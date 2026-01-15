function X = Secant(distribution_type, parameters, delta, n, x_min, x_max)
    X = zeros(1, n);
    F = @(x) ContinuousCDF(x, distribution_type, parameters);
    if(n < 1 || delta <= 0)
       error("Wrong parameters"); 
    end
    umin = F(x_min); umax = F(x_max);
    for i = 1:n
       U = ULEcuyerRNG(1) * (umax - umin) + umin;
       a = x_min; b = x_max;
       x = a + (b - a) * (U - F(a)) / (F(b) - F(a));
       while(b - a > delta && abs(U - F(x)) > delta)
          if(U < F(x))
              b = x;
          else
              a = x;
          end
        x = a + (b - a) * (U - F(a)) / (F(b) - F(a));
       end
       X(i) = x;
    end
end