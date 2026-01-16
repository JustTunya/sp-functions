function [X, Y] = CompBoundInv(distribution_type, intervals, n, M)
    ax = intervals(1, 1); bx = intervals(1, 2);
    ay = intervals(2, 1); by = intervals(2, 2);

    f = @(x, y) ContinuousPDF_2D(x, y, distribution_type, []);

    X = zeros(1, n); Y = zeros(1, n);

    for i = 1:n
        while true
            U = ULEcuyerRNG(1);
            x = ax + (bx - ax) * ULEcuyerRNG(1);
            y = ay + (by - ay) * ULEcuyerRNG(1);

            if (U * M <= f(x, y))
                X(i) = x;
                Y(i) = y;
                break;
            end
        end
    end
end