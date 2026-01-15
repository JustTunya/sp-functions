function [X1, X2] = MarsagliaStandNorm(n)
    X1 = zeros(1, n); X2 = zeros(1, n);

    for i = 1:n
        while (true)
            U1 = ULEcuyerRNG(1); U2 = ULEcuyerRNG(1);

            Z1 = 2 * U1 - 1; Z2 = 2 * U2 - 1;
            S = Z1^2 + Z2^2;

            if (S > 0 && S <= 1)
                T = sqrt(-2 * log(S) / S);
                X1(i) = T * Z1;
                X2(i) = T * Z2;
                break;
            end
        end
    end
end