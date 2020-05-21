function [tOut, y] = bdf1si(odefun, tOut, yn, options)
    if nargin < 4
        options = odeset();
    end
    options = ensureSolverSet(options);

    n = length(tOut);
    y = zeros(length(yn), n);

    [t, saveIndices] = timepointsWithMaxStep(tOut, options);

    y(:, 1) = yn;
    j = 2;

    yCoeff = [1, -1]';
    implicitCoeff = -1;
    explicitCoeff = -1;

    for i = 2:length(t)
        dt = t(i) - t(i-1);

        explicitF = yn * yCoeff(2) / dt + ...
            odefun.explicit(t(i-1), yn) * explicitCoeff;

        yn = options.optimmethod( ...
            @(h) fun(odefun.implicit, t(i), h, dt, explicitF, options), ...
            yn, ...
            options.optimoptions);

        if i == saveIndices(j)
            y(:,j) = yn;
            j = j + 1;
        end

        if any(isnan(yn))
            fprintf('Nan`s in solution\n')
            break;
        end
    end

    function [F, J] = fun(implicitOdefun, t, h, dt, explicitF, options)
        f = implicitOdefun(t, h);
        F = h * yCoeff(1) / dt + f * implicitCoeff + explicitF;

        if nargout == 2
            jac = options.Jacobian(t, h);
            J = speye(length(f)) * yCoeff(1) / dt + jac * implicitCoeff;
        end
    end

    [tOut, y] = functionOutputParser(tOut, y, nargout);
end

