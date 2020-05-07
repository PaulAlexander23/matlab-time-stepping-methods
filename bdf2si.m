function [t, y] = bdf2si(odefun, t, y0, options)
    if nargin < 4
        optimopts = optimoptions('fsolve', 'Display', 'off');
        options = struct('optimmethod', @fsolve, 'optimoptions', optimopts);
    end

    n = length(t);
    y = zeros(length(y0), n);

    windup = bdf1si(odefun, t(1:2), y0, options);
    y(:, 1:2) = windup.y';

    yCoeff = [1, -4/3, 1/3]';
    explicitCoeff = -[4/3, -2/3]';
    implicitCoeff = -2/3;

    for i = 3:n
        dt = t(i) - t(i-1);

        explicitF = y(:, i-1:-1:i-2) * yCoeff(2:3) + ...
            dt * [odefun.explicit(t(i-1), y(:, i-1)), ...
            odefun.explicit(t(i-2), y(:, i-2))] * explicitCoeff;

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
