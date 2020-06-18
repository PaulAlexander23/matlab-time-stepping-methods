function [tOut, y] = bdf4si(odefun, tOut, y0, options)
    if nargin < 4
        options = odeset();
    end
    options = ensureSolverSet(options);

    hasEvents = ~isempty(options.Events);

    n = length(tOut);
    y = zeros(size(y0, 1), n);
    quit = false;
    if hasEvents, value = 1; end

    [t, saveIndices] = timepointsWithMaxStep(tOut, options);
    validateTimeStepsEqual(t);

    yCoeff = [25, -48, 36, -16, 3]'/25;
    explicitCoeff = -[48, -72, 48, -12]'/25;
    implicitCoeff = -12/25;

    if size(y0, 2) == 1
        windup = bdf3si(odefun, t(1:4), y0, options);
        y0 = windup.y';
    end

    i = 0;
    k = 1;
    while ~quit && i < 4
        i = i + 1;

        if i >= 4, ynm3 = ynm2; end
        if i >= 3, ynm2 = ynm1; end
        if i >= 2, ynm1 = yn; end

        yn = y0(:,i);

        if i == saveIndices(k)
            y(:,k) = yn;
            k = k + 1;
        end

        if hasEvents
            [quit, value, ie, xe, ye] = handleEvents(options.Events, t(i), yn, value);
        end
    end


    while ~quit
        i = i + 1;

        dt = t(i) - t(i-1);

        explicitF = [yn, ynm1, ynm2, ynm3] * yCoeff(2:end) + ...
            dt * [odefun.explicit(t(i-1), yn), ...
            odefun.explicit(t(i-2), ynm1), ...
            odefun.explicit(t(i-3), ynm2), ...
            odefun.explicit(t(i-4), ynm3)] * explicitCoeff;

        ynm3 = ynm2;
        ynm2 = ynm1;
        ynm1 = yn;

        yn = options.optimmethod( ...
            @(h) fun(odefun.implicit, t(i), h, dt, explicitF, options), ...
            yn, ...
            options.optimoptions);

        if i == saveIndices(k)
            y(:,k) = yn;
            k = k + 1;
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
        F = h * yCoeff(1) + dt * f * implicitCoeff + explicitF;

        if nargout == 2
            j = options.Jacobian(t, h);
            J = speye(length(f)) * yCoeff(1) + dt * j * implicitCoeff;
        end

    end

    if ~hasEvents
        [tOut, y] = functionOutputParser(tOut, y, nargout);
    else
        [tOut, y] = functionOutputParserEvents(tOut, y, ie, xe, ye, nargout);
    end
end
