function [tOut, y] = ab2(odefun, tOut, yn, options)
    if nargin < 4, options = odeset(); end

    n = length(tOut);
    y = zeros(length(yn),n);

    [t, saveIndices] = timepointsWithMaxStep(tOut, options);
    validateTimeStepsEqual(t);
    
    y(:,1) = yn;
    j = 2;

    
    for i = 2:length(t)
        if i == 2
            ynm1 = yn;
            yn = yn + (t(2)-t(1)) * odefun(t(1), yn);
        else
            ytemp = yn;
            yn = yn + (t(i)-t(i-1)) * (3/2 * odefun(t(i-1), yn) - 1/2 * odefun(t(i-2), ynm1));
            ynm1 = ytemp;
        end

        if i == saveIndices(j)
            y(:,j) = yn;
            j = j + 1;
        end
    end

    [tOut, y] = functionOutputParser(tOut, y, nargout);
end
