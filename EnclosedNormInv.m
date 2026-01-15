function[X, Sample, Density] = EnclosedNormInv(n)
    if n < 1
        error("The parameter is wrong!");
    end

    alpha = 1 / sqrt(exp(1));
    beta = 1/2;
    gamma = sqrt(2);

    X = zeros(1, n);
    n_density = 0; n_sample = 0;

    for i = 1:n
        while true
            n_density = n_density + 1;

            U = ULEcuyerRNG(1);
            V = ULEcuyerRNG(1);
            Y = tan(pi * V);
            S = beta * (Y^2);
            W = (alpha * U) / (beta + S);
    
            if (abs(Y) > gamma)
                L = false;
            else
                L = (W <= 1 - S);
            end

            if ~L
                n_sample = n_sample + 1;
                L = (W <= exp(-S));
            end
            
            if L
                break;
            end
        end

        X(i) = Y;
    end

    Density = n_density / n;
    Sample = n_sample / n;
end