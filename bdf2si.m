function [t, y] = bdf2si(explicitOdefun, implicitOdefun, t, y0, options)
    if nargin < 4
        options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
            optimoptions('fsolve', 'Display', 'off')));
    end

    n = length(t);
    y = zeros(length(y0), n);

    windup = bdf1si(explicitOdefun, implicitOdefun, t(1:2), y0, options);
    y(:, 1:2) = windup.y';

    yCoeff = [1, -4/3, 1/3]';
    explicitCoeff = -[4/3, -2/3]';
    implicitCoeff = -2/3;

    for i = 3:n
        dt = t(i) - t(i-1);

        explicitF = y(:, i-1:-1:i-2) * yCoeff(2:3) + ...
            dt * [explicitOdefun(t(i-1), y(:, i-1)), ...
            explicitOdefun(t(i-2), y(:, i-2))] * explicitCoeff;

        y(:, i) = options.optimmethod(@(h) fun(implicitOdefun, t(i), h, dt, explicitF), ...
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
