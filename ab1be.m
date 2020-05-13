function [tOut, y] = ab1be(odefun, tOut, yn, options)
    if nargin < 4
        options = odeset();
    end
    options = ensureSolverSet(options);

    n = length(tOut);
    y = zeros(length(yn),n);
    
    [t, saveIndices] = timepointsWithMaxStep(tOut, options);

    y(:,1) = yn;
    j = 2;

    explicitCoeff = 1;
    implicitCoeff = 1;
    
    for i = 2:length(t)
        ynm1 = yn;
        dt = t(i) - t(i - 1);
        yn = ynm1 + ...
            dt * odefun.explicit(t(i-1), ynm1) * explicitCoeff;
        yn = options.optimmethod( ...
            @(x) x - yn - dt * odefun.implicit(t(i), x) * implicitCoeff,...
            ynm1, ...
            options.optimoptions);

        if i == saveIndices(j)
            y(:,j) = yn;
            j = j + 1;
        end
    end

    [tOut, y] = functionOutputParser(tOut, y, nargout);
end
