function [t, y] = bdf3si(odefun, t, y0, options)
    if nargin < 4
        options = odeset();
    end

    windup = bdf2si(odefun, t(1:3), y0, options);

    [t, y] = timesStepperInterface(odefun, t, y0, options, scheme, windup);

    function y = scheme(i, t, y, options)
        yCoeff = [11, -18, 9, -2]'/11;
        explicitCoeff = -[18, -18, 6]'/11;
        implicitCoeff = -6/11;

        dt = t(i) - t(i-1);

        explicitF = y(:, i-1:-1:i-3) * yCoeff(2:end) + ...
            dt * [odefun.explicit(t(i-1), y(:, i-1)), ...
            odefun.explicit(t(i-2), y(:, i-2)), ...
            odefun.explicit(t(i-3), y(:, i-3))] * explicitCoeff;

        y(:, i) = options.optimmethod( ...
            @(h) fun(odefun.implicit, t(i), h, dt, explicitF, options), ...
            y(:, i - 1), ...
            options.optimoptions);

        function [F, J] = fun(implicitOdefun, t, h, dt, explicitF, options)
            f = implicitOdefun(t, h);
            F = h * yCoeff(1) + dt * f * implicitCoeff + explicitF;

            if nargout == 2
                j = options.Jacobian(t, h);
                J = speye(length(f)) * yCoeff(1) + dt * j * implicitCoeff;
            end
        end
    end
end
