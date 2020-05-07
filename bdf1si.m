function [t, y] = bdf1si2(odefun, t, y0, options)
    if nargin < 4
        options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
            optimoptions('fsolve', 'Display', 'off')));
    end

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

