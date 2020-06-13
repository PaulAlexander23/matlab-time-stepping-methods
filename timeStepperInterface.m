function timeStepperInterface(odefun, t, y0, options, scheme, windup)
    if nargin < 4
        options = odeset();
    end
    options = ensureSolverSet(options);

    n = length(t);
    y = zeros(size(y0, 1), n);

    if size(y0, 2) == 1
        y(:, 1:3) = windup();
    else
        y(:, 1:3) = y0;
    end

    for i = 4:n
        y(:, i) = scheme(t, y);

        if any(isnan(y(:, i)))
            fprintf('Nan`s in solution\n')
            break;
        end
    end

    [t, y] = functionOutputParser(t, y, nargout);
end
