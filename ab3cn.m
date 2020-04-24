function [t, y] = ab3cn(explicitOdefun, implicitOdefun, t, y0, options)
    if nargin < 4
        options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
            optimoptions('fsolve', 'Display', 'off')));
    end

    n = length(t);
    y = zeros(length(y0),n);
    
    windup = ab2be(explicitOdefun, implicitOdefun, t(1:3), y0, options);
    y(:, 1:3) = windup.y';
    
    explicitCoeff = [23/12,-4/3,5/12]';
    implicitCoeff = [3/2,-1/2]';
    
    for i = 4:n
        dt = t(i) - t(i - 1);
        y(:, i) = y(:, i-1) + ...
            dt * [explicitOdefun(t(i-1), y(:, i-1)), ...
                explicitOdefun(t(i-2), y(:, i-2)), ...
                explicitOdefun(t(i-3), y(:, i-3))] *  explicitCoeff;
        y(:, i) = options.optimmethod(@(x) x - y(:, i) - ...
            dt * [implicitOdefun(t(i), x), ...
                implicitOdefun(t(i-1), y(:, i-1))] * implicitCoeff, ...
            y(:, i-1));
    end

    [t, y] = functionOutputParser(t, y, nargout);
end
