function [tOut, y] = bdf1si(odefun, tOut, yn, options)
    if nargin < 4
        options = odeset();
    end
    options = ensureSolverSet(options);

    hasEvents = ~isempty(options.Events);

    n = length(tOut);
    y = zeros(length(yn), n);
    quit = false;
    if hasEvents, value = 1; end

    [t, saveIndices] = timepointsWithMaxStep(tOut, options);

    yCoeff = [1, -1]';
    implicitCoeff = -1;
    explicitCoeff = -1;

    i = 1;
    y(:, i) = yn;
    j = 2;
    if hasEvents
        [~, value] = handleEvents(options.Events, t(i), yn, value);
    end

    while ~quit 
        i = i + 1;

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

        if hasEvents
            [quit, value, ie, xe, ye] = handleEvents(options.Events, t(i), yn, value);
        end

        if i >= length(t)
            quit = true;
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

    if ~hasEvents
        [tOut, y] = functionOutputParser(tOut, y, nargout);
    else
        [tOut, y] = functionOutputParserEvents(tOut, y, ie, xe, ye, nargout);
    end
end

