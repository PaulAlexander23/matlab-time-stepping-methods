function [tOut, y] = ab2be(odefun, tOut, yn, options)
    if nargin < 4
        options = odeset();
    end
    options = ensureSolverSet(options);

    n = length(tOut);
    y = zeros(length(yn),n);

    [t, saveIndices] = timepointsWithMaxStep(tOut, options);
    validateTimeStepsEqual(t);

    y(:,1) = yn;
    j = 2;

    explicitCoeff = [3/2,-1/2]';
    implicitCoeff = 1;

    for i = 2:length(t)
        dt = t(i) - t(i - 1);
        if i == 2
            ynm1 = yn;
            yn = ynm1 + ...
                dt * odefun.explicit(t(i-1), ynm1);
            yn = options.optimmethod( ...
                @(x) x - yn - dt * odefun.implicit(t(i), x),...
                ynm1, ...
                options.optimoptions);
        else
            ynm2 = ynm1;
            ynm1 = yn;
            yn = ynm1 + ...
                dt * [odefun.explicit(t(i-1), ynm1), ...
                odefun.explicit(t(i-2), ynm2)] * explicitCoeff;
            yn = options.optimmethod( ...
                @(x) x - yn - dt * odefun.implicit(t(i), x) * implicitCoeff,...
                ynm1, ...
                options.optimoptions);
        end

        if i == saveIndices(j)
            y(:,j) = yn;
            j = j + 1;
        end
    end

    [tOut, y] = functionOutputParser(tOut, y, nargout);
end
