function [tOut, y] = bdf2si(odefun, tOut, yn, options)
    if nargin < 4
        options = odeset();
    end
    options = ensureSolverSet(options);

    n = length(tOut);
    y = zeros(length(yn), n);

    [t, saveIndices] = timepointsWithMaxStep(tOut, options);
    validateTimeStepsEqual(t);

    y(:, 1) = yn;
    j = 2;

    yCoeff = [1, -4/3, 1/3]';
    explicitCoeff = -[4/3, -2/3]';
    implicitCoeff = -2/3;

    for i = 2:length(t)
        dt = t(i) - t(i-1);
        if i == 2
            explicitF = - yn - ...
                dt * odefun.explicit(t(i-1), yn);
            ynm1 = yn;
            yn = options.optimmethod( ...
                @(h) funBDF1(odefun.implicit, t(i), h, dt, explicitF, options), ...
                yn, ...
                options.optimoptions);

        else
            explicitF = [yn, ynm1] * yCoeff(2:3) + ...
                dt * [odefun.explicit(t(i-1), yn), ...
                odefun.explicit(t(i-2), ynm1)] * explicitCoeff;
            ynm1 = yn;
            yn = options.optimmethod( ...
                @(h) funBDF2(odefun.implicit, t(i), h, dt, explicitF, options), ...
                ynm1, ...
                options.optimoptions);
        end

        if i == saveIndices(j)
            y(:,j) = yn;
            j = j + 1;
        end

        if any(isnan(yn))
            fprintf('Nan`s in solution\n')
            break;
        end
    end

    function [F, J] = funBDF1(implicitOdefun, t, h, dt, explicitF, options)
        f = implicitOdefun(t, h);
        F = h - dt * f + explicitF;

        if nargout == 2
            jac = options.Jacobian(t, h);
            J = speye(length(f)) - dt * jac;
        end
    end
    
    function [F, J] = funBDF2(implicitOdefun, t, h, dt, explicitF, options)
        f = implicitOdefun(t, h);
        F = h * yCoeff(1) + dt * f * implicitCoeff + explicitF;
        if nargout == 2
            jac = options.Jacobian(t, h);
            J = speye(length(f)) * yCoeff(1) + dt * jac * implicitCoeff;
        end

    end

    [tOut, y] = functionOutputParser(tOut, y, nargout);
end
