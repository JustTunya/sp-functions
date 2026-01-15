function[X, FEvals, TIters] = EnclosedCustom(n)
    if n < 1
        error("The parameter is wrong!");
    end

    a = 0; b = pi / 3;
    alpha = 2 / pi;
    g = 3 / pi;
    h1 = 5 / (2 * pi);
    h2 = 7 / (2 * pi);
    c = h2 / g;

    X = zeros(1);
    n_t_iters = 0;
    n_f_evals = 0;

    for i  = 1:n
        while true
            n_t_iters = n_t_iters + 1;

            U = ULEcuyerRNG;
            [Yvec, ~] = URealRNG(1, "ULEcuyer", a, b, 1);
            Y = Yvec(1);
            W = U * c * g;

            if (W <= h1)
                X(i) = Y;
                break;
            end

            n_f_evals = n_f_evals + 1;
            fY = 1/2 * alpha * (sin(3 * Y)^2 + 5 / 2);

            if W <= fY
                X(i) = Y;
                break
            end
        end
    end

    FEvals = n_f_evals / n;
    TIters = n_t_iters / n;
end