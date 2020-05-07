function [t, y] = bdf1si(odefun, t, y0, options)
    if nargin < 4
        options = odeset();
    end
    options = ensureSolverSet(options);

    n = length(t);
    y = zeros(length(y0), n);

    y(:, 1) = y0;

    yCoeff = [1, -1]';
    implicitCoeff = -1;
    explicitCoeff = -1;

    for i = 2:n
        dt = t(i) - t(i-1);

        explicitF = y(:, i-1) * yCoeff(2) + ...
            dt * odefun.explicit(t(i-1), y(:, i-1)) * explicitCoeff;

        y(:, i) = options.optimmethod( ...
            @(h) fun(odefun.implicit, t(i), h, dt, explicitF, options), ...
            y(:, i - 1), ...
            options.optimoptions);

        if any(isnan(y(:, i)))
            fprintf('Nan`s in solution\n')
            break;
        end
    end

    function [F, J] = fun(implicitOdefun, t, h, dt, explicitF, options)
        f = implicitOdefun(t, h);
        F = h * yCoeff(1) + dt * f * implicitCoeff + explicitF;

        if nargout == 2
            j = options.Jacobian(t, h);
            J = speye(length(f)) * yCoeff(1) + dt * j * implicitCoeff;
        end
    end

    [t, y] = functionOutputParser(t, y, nargout);
end

