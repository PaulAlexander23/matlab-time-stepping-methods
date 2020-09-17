function [tOut, y] = bdf4si(odefun, tOut, y0, options)
    if nargin < 4
        options = odeset();
    end
    options = ensureSolverSet(options);

    hasEvents = ~isempty(options.Events);

    if isfield(options, "MaxStep")
        if ~isempty(options.MaxStep)
            dt = options.MaxStep;
        else
            dt = tOut(2);
        end
    else
        dt = tOut(2);
    end

    n = length(tOut);
    y = zeros(size(y0, 1), n);
    quit = false;
    if hasEvents, value = 1; end

    yCoeff = [25, -48, 36, -16, 3]'/25;
    explicitCoeff = -[48, -72, 48, -12]'/25;
    implicitCoeff = -12/25;

    if size(y0, 2) == 1
        windup = bdf3si(odefun, dt * (0:3)', y0, options);
        y0 = windup.y';
    end

    i = 1;
    k = 1;
    t = 0;
    while ~quit
        if i < 5
            if i >= 4, ynm3 = ynm2; end
            if i >= 3, ynm2 = ynm1; end
            if i >= 2, ynm1 = yn; end

            yn = y0(:,i);
        else
            explicitF = [yn, ynm1, ynm2, ynm3] * yCoeff(2:end) + ...
                dt * [odefun.explicit(t-dt, yn), ...
                odefun.explicit(t - 2*dt, ynm1), ...
                odefun.explicit(t - 3*dt, ynm2), ...
                odefun.explicit(t - 4*dt, ynm3)] * explicitCoeff;

            ynm3 = ynm2;
            ynm2 = ynm1;
            ynm1 = yn;

            yn = options.optimmethod( ...
                @(h) fun(odefun.implicit, t, h, dt, explicitF, options), ...
                yn, ...
                options.optimoptions);
        end

        if abs(t - tOut(k)) < 1e-13
            y(:,k) = yn;
            k = k + 1;
        end

        if hasEvents
            [quit, value, ie, xe, ye] = handleEvents(options.Events, t, yn, value);
            if i == 1, quit = false; end
        end

        if abs(t - tOut(end)) < 1e-13
            quit = true;
        end

        i = i + 1;
        t = t + dt;
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
