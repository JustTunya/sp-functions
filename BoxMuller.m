function [X1, X2] = BoxMuller(n)
    U = URealRNG(1, "ULEcuyer", 0, 1, [2, n]);
    R = sqrt(-2 * log(U(1,:)));
    Theta = 2 * pi * U(2,:);

    X1 = R .* cos(Theta);
    X2 = R .* sin(Theta);
end