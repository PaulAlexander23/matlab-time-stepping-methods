function y = ab1be(explicitOdefun, implicitOdefun, t, y0, options)
    if nargin < 4
        options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
            optimoptions('fsolve', 'Display', 'off')));
    end

    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;

    explicitCoeff = 1;
    implicitCoeff = 1;
    
    for i = 2:n
        dt = t(i) - t(i - 1);
        y(:, i) = y(:, i-1) + ...
            dt * explicitOdefun(t(i-1), y(:, i-1)) * explicitCoeff;
        y(:, i) = options.optimmethod(@(x) x - y(:, i) - ...
            dt * implicitOdefun(t(i), x) * implicitCoeff, y(:, i-1));
    end
end