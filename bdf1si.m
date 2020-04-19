function [t, y] = bdf1si(explicitOdefun, implicitOdefun, t, y0, options)
    if nargin < 4
        options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
            optimoptions('fsolve', 'Display', 'off')));
    end

    n = length(t);
    y = zeros(length(y0), n);

    y(:, 1) = y0;

    explicitCoeff = 1;
    implicitCoeff = [1, -1]';

    for i = 2:n
        explicitF = explicitOdefun(t(i), y(:, i-1)) * explicitCoeff;
        y(:, i) = options.optimmethod(@(h) fun(h, implicitOdefun, i, t, y, explicitF), ...
            y(:, i - 1));

        if any(isnan(y(:, i)))
            fprintf('Nan`s in solution\n')
            break;
        end
    end

    function [F, J] = fun(h, implicitOdefun, i, t, y, expf)
        dt = t(i) - t(i-1);

        if nargout == 1
            f = implicitOdefun(t(i), h);
        elseif nargout == 2
            [f, j] = implicitOdefun(t(i), h);
            J = implicitCoeff(1) * speye(length(f)) / dt - j;
        end
        
        F = ([h, y(:, i-1)] * implicitCoeff) / dt - f;
        F = F - expf;
    end

    [t, y] = functionOutputParser(t, y, nargout);
end
