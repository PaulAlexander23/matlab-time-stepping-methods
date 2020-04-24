function [t, y] = bdf3si(explicitOdefun, implicitOdefun, t, y0, options)
    if nargin < 5
        options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
            optimoptions('fsolve', 'Display', 'off')));
    end

    n = length(t);
    y = zeros(length(y0), n);

    [~, y0] = bdf2si(explicitOdefun, implicitOdefun, t(1:3), y0, options);
    y(:, 1:3) = y0';

    explicitCoeff = [2, -1]';
    implicitCoeff = [11/6, -3, 3/2, -1/3]';

    for i = 4:n
        explicitF = [explicitOdefun(t(i-1), y(:, i-1)), ...
            explicitOdefun(t(i-2), y(:, i-2))] * explicitCoeff;
        y(:, i) = options.optimmethod(@(h) fun(h, implicitOdefun, i, t, y, explicitF), ...
            y(:, i - 1));

        if any(isnan(y(:, i)))
            fprintf('Nan`s in solution\n')
            break;
        end
    end

    function [F, J] = fun(h, implicitOdefun, i, t, y, explicitF)
        dt = t(i) - t(i-1);

        if nargout == 1
            f = implicitOdefun(t(i), h);
        elseif nargout == 2
            [f, j] = implicitOdefun(t(i), h);
            J = implicitCoeff(1) * speye(length(f)) / dt -  j;
        end

        F = ([h, y(:, i-1:-1:i-3)] * implicitCoeff) / dt - f;
        F = F - explicitF;
    end

    [t, y] = functionOutputParser(t, y, nargout);
end
