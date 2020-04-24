function [t, y] = ab2be(explicitOdefun, implicitOdefun, t, y0, options)
    if nargin < 4
        options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
            optimoptions('fsolve', 'Display', 'off')));
    end

    n = length(t);
    y = zeros(length(y0),n);

    windup = ab1be(explicitOdefun, implicitOdefun, t(1:2), y0, options);
    y(:,1:2) = windup.y';

    explicitCoeff = [3/2,-1/2]';
    implicitCoeff = 1;

    for i = 3:n
        dt = t(i) - t(i - 1);
        y(:, i) = y(:, i-1) + ...
            dt * [explicitOdefun(t(i-1), y(:, i-1)), ...
            explicitOdefun(t(i-2), y(:, i-2))] * explicitCoeff;
        y(:, i) = options.optimmethod(@(x) x - y(:, i) - ...
            dt * implicitOdefun(t(i), x) * implicitCoeff, y(:, i-1));
    end

    [t, y] = functionOutputParser(t, y, nargout);
end
