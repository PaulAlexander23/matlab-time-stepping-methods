function [t, y] = bdf4si(odefun, t, y0, options)
    if nargin < 5
        options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
            optimoptions('fsolve', 'Display', 'off')));
    end

    n = length(t);
    y = zeros(length(y0), n);

    windup = bdf3si(odefun, t(1:3), y0, options);
    y(:, 1:3) = windup.y';

    yCoeff = [25, -48, 36, -16, 3]'/25;
    explicitCoeff = -[48, -72, 48, -12]'/25;
    implicitCoeff = -12/25;

    for i = 5:n
        dt = t(i) - t(i-1);

        explicitF = y(:, i-1:-1:i-4) * yCoeff(2:end) + ...
            dt * [odefun.explicit(t(i-1), y(:, i-1)), ...
            odefun.explicit(t(i-2), y(:, i-2)), ...
            odefun.explicit(t(i-3), y(:, i-3)), ...
            odefun.explicit(t(i-4), y(:, i-4))] * explicitCoeff;

        y(:, i) = options.optimmethod(@(h) fun(odefun.implicit, t(i), h, dt, explicitF), ...
            y(:, i - 1));

        if any(isnan(y(:, i)))
            fprintf('Nan`s in solution\n')
            break;
        end
    end

    function [F, J] = fun(implicitOdefun, t, h, dt, explicitF)
        if nargout == 1
            f = implicitOdefun(t, h);
        elseif nargout == 2
            [f, j] = implicitOdefun(t, h);
            J = speye(length(f)) * yCoeff(1) + dt * j * implicitCoeff;
        end

        F = h * yCoeff(1) + dt * f * implicitCoeff + explicitF;
    end

    [t, y] = functionOutputParser(t, y, nargout);
end


