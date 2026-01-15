function X = NewtonRaphson(distribution_type, parameters, urng, limits, eps, n, initial_value)
    U = genU(n, urng);
    U = min(max(U, 1e-12), 1-1e-12);

    X = zeros(1,n);
    switch string(distribution_type)
        case "gamma"
            a = parameters(1); b = parameters(2);
            if ~(a>0 && b>0), error("gamma: a,b > 0 kell"); end

            if nargin < 7 || isempty(initial_value) || ~isfinite(initial_value)
                if a > 1
                    x0 = (a-1)*b;
                else
                    x0 = a*b;
                end
            else
                x0 = initial_value;
            end
            x0 = min(max(x0, limits(1)), limits(2));

            F = @(x) ContinuousCDF(x, distribution_type, parameters);
            f = @(x) ContinuousPDF(x, distribution_type, parameters);

            maxit = 200;
            for i = 1:n
                u = U(i);
                x = x0;
                for it = 1:maxit
                    Fx = F(x);
                    if abs(Fx - u) <= eps, break; end
                    d = f(x);
                    if d <= 0 || ~isfinite(d)
                        if Fx < u
                            x = 0.5*(x + limits(2));
                        else
                            x = 0.5*(x + limits(1));
                        end
                        continue
                    end
                    xNew = x - (Fx - u)/d;
                    if xNew <= limits(1) || xNew >= limits(2) || ~isfinite(xNew)
                        if Fx < u
                            xNew = 0.5*(x + limits(2));
                        else
                            xNew = 0.5*(x + limits(1));
                        end
                    end
                    if abs(xNew - x) <= max(1e-12, eps*max(1,x)), x = xNew; break; end
                    x = xNew;
                end
                X(i) = min(max(x, limits(1)), limits(2));
            end
        case "L12"
            theta = parameters(1);
            if ~(theta > 0)
                error("L12: theta must be strictly positive!");
            end

            a = limits(1); b = limits(2);
            if ~isfinite(a) || ~isfinite(b) || a >= b
                error("L12: invalid limits!");
            end

            a = max(a, 0);
            b = min(b, 2);
            if a >= b
                error("L12: limits must intersect [0,2] with positive length.");
            end

            maxit = 200;

            for i = 1:n
                u = U(i);

                if nargin >= 7 && ~isempty(initial_value) && isfinite(initial_value)
                    x = initial_value;
                else
                    x = 2 * (1 - (1 - u)^(1/(3*theta)));
                end
                x = min(max(x, a), b);

                lo = a; hi = b;

                for it = 1:maxit
                    if x <= 0
                        Fx = 0;
                        d  = (3*theta/2);
                    elseif x >= 2
                        Fx = 1;
                        d  = 0;
                    else
                        t = 1 - x/2;
                        Fx = 1 - t^(3*theta);
                        d  = (3*theta/2) * t^(3*theta - 1);
                    end

                    err = Fx - u;
                    if abs(err) <= eps
                        break;
                    end

                    if Fx < u
                        lo = max(lo, x);
                    else
                        hi = min(hi, x);
                    end

                    if d <= 0 || ~isfinite(d)
                        xNew = 0.5 * (lo + hi);
                    else
                        xNew = x - err / d;

                        if ~isfinite(xNew) || xNew <= lo || xNew >= hi
                            xNew = 0.5 * (lo + hi);
                        end
                    end

                    if abs(xNew - x) <= max(1e-12, eps*max(1,abs(x)))
                        x = xNew;
                        break;
                    end
                    x = xNew;
                end

                X(i) = min(max(x, a), b);
            end
        otherwise
            error("Nem támogatott eloszlás: %s", distribution_type);
    end
end

function U = genU(m, urng)
    u = lower(string(urng));
    switch u
        case {"ulecuyer","lecuyer","1"}
            U = ULEcuyerRNG(m);
        case {"mersennetwister","umersennetwister","2"}
            U = UMersenneTwisterRNG(m);
        otherwise
            U = rand(1,m);
    end
end