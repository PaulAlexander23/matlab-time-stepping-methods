function [tOut, y] = bdf2si(odefun, tOut, y0, options)
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

    yCoeff = [1, -4/3, 1/3]';
    explicitCoeff = -[4/3, -2/3]';
    implicitCoeff = -2/3;

    if size(y0, 2) == 1
        windup = bdf1si(odefun, t(1:2), y0, options);
        y0 = windup.y';
    end

    i = 0;
    k = 1;
    while ~quit && i < 2
        i = i + 1;

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
        if i == 2
            explicitF = - yn / dt - ...
                odefun.explicit(t(i-1), yn);
            ynm1 = yn;
            yn = options.optimmethod( ...
                @(h) funBDF1(odefun.implicit, t(i), h, dt, explicitF, options), ...
                yn, ...
                options.optimoptions);

        else
            explicitF = [yn, ynm1] * yCoeff(2:3) / dt + ...
                [odefun.explicit(t(i-1), yn), ...
                odefun.explicit(t(i-2), ynm1)] * explicitCoeff;
            ynm1 = yn;
            yn = options.optimmethod( ...
                @(h) funBDF2(odefun.implicit, t(i), h, dt, explicitF, options), ...
                ynm1, ...
                options.optimoptions);
        end

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

    function [F, J] = funBDF1(implicitOdefun, t, h, dt, explicitF, options)
        f = implicitOdefun(t, h);
        F = h / dt - f + explicitF;

        if nargout == 2
            jac = options.Jacobian(t, h);
            J = speye(length(f)) / dt - jac;
        end
    end
    
    function [F, J] = funBDF2(implicitOdefun, t, h, dt, explicitF, options)
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
