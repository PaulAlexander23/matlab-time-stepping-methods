function [t, y] = ab1be(odefun, t, y0, options)
    if nargin < 4
        optimopts = optimoptions('fsolve', 'Display', 'off');
        options = struct('optimmethod', @fsolve, 'optimoptions', optimopts);
    end

    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;

    explicitCoeff = 1;
    implicitCoeff = 1;
    
    for i = 2:n
        dt = t(i) - t(i - 1);
        y(:, i) = y(:, i-1) + ...
            dt * odefun.explicit(t(i-1), y(:, i-1)) * explicitCoeff;
        y(:, i) = options.optimmethod( ...
            @(x) x - y(:, i) - dt * odefun.implicit(t(i), x) * implicitCoeff,...
            y(:, i-1), ...
            options.optimoptions);
    end

    [t, y] = functionOutputParser(t, y, nargout);
end
