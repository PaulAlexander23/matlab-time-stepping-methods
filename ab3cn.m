function [t, y] = ab3cn(odefun, t, y0, options)
    if nargin < 4
        optimopts = optimoptions('fsolve', 'Display', 'off');
        options = struct('optimmethod', @fsolve, 'optimoptions', optimopts);
    end

    n = length(t);
    y = zeros(length(y0),n);
    
    windup = ab2be(odefun, t(1:3), y0, options);
    y(:, 1:3) = windup.y';
    
    explicitCoeff = [23/12,-4/3,5/12]';
    implicitCoeff = [3/2,-1/2]';
    
    for i = 4:n
        dt = t(i) - t(i - 1);
        y(:, i) = y(:, i-1) + ...
            dt * [odefun.explicit(t(i-1), y(:, i-1)), ...
                odefun.explicit(t(i-2), y(:, i-2)), ...
                odefun.explicit(t(i-3), y(:, i-3))] *  explicitCoeff;
        y(:, i) = options.optimmethod( ...
            @(x) x - y(:, i) - ...
            dt * [odefun.implicit(t(i), x), ...
                odefun.implicit(t(i-1), y(:, i-1))] * implicitCoeff, ...
            y(:, i-1), ...
            options.optimoptions);
    end

    [t, y] = functionOutputParser(t, y, nargout);
end
