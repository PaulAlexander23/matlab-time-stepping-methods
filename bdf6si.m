function [t, y] = bdf6si(odefun, t, y0, options)
    if nargin < 4
        options = odeset();
    end
    options = ensureSolverSet(options);

    n = length(t);
    y = zeros(size(y0, 1), n);

    if size(y0, 2) == 1
        windup = bdf5si(odefun, t(1:6), y0, options);
        y(:, 1:6) = windup.y';
    else
        y(:, 1:6) = y0;
    end

    yCoeff = [147/60, -6, 15/2, -20/3, 15/4, -6/5, 1/6]';
    explicitCoeff = -[6, -15, 20, -15, 6, -1]';
    implicitCoeff = -1;

    for i = 7:n
        dt = t(i) - t(i-1);

        explicitF = y(:, i-1:-1:i-6) * yCoeff(2:end) + ...
            dt * [odefun.explicit(t(i-1), y(:, i-1)), ...
            odefun.explicit(t(i-2), y(:, i-2)), ...
            odefun.explicit(t(i-3), y(:, i-3)), ...
            odefun.explicit(t(i-4), y(:, i-4)), ...
            odefun.explicit(t(i-5), y(:, i-5)), ...
            odefun.explicit(t(i-6), y(:, i-6))] * explicitCoeff;

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
