function [tOut, y] = bdf6si(odefun, tOut, y0, options)
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

    yCoeff = [147/60, -6, 15/2, -20/3, 15/4, -6/5, 1/6]';
    explicitCoeff = -[6, -15, 20, -15, 6, -1]';
    implicitCoeff = -1;

    if size(y0, 2) == 1
        windup = bdf5si(odefun, t(1:6), y0, options);
        y0 = windup.y';
    end

    i = 0;
    k = 1;
    while ~quit && i < 6
        i = i + 1;

        if i >= 6, ynm5 = ynm4; end
        if i >= 5, ynm4 = ynm3; end
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
            if i == 1, quit = false; end
        end
    end


    while ~quit
        i = i + 1;

        dt = t(i) - t(i-1);

        explicitF = [yn, ynm1, ynm2, ynm3, ynm4, ynm5] * yCoeff(2:end) + ...
            dt * [odefun.explicit(t(i-1), yn), ...
            odefun.explicit(t(i-2), ynm1), ...
            odefun.explicit(t(i-3), ynm2), ...
            odefun.explicit(t(i-4), ynm3), ...
            odefun.explicit(t(i-5), ynm4), ...
            odefun.explicit(t(i-6), ynm5)] * explicitCoeff;

        ynm5 = ynm4;
        ynm4 = ynm3;
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
