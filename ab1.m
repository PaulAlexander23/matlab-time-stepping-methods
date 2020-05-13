function [tOut, y] = ab1(odefun, tOut, y0, options)
    if nargin < 4, options = odeset(); end

    n = length(tOut);
    y = zeros(length(y0),n);
    
    [t, saveIndices] = timepointsWithMaxStep(tOut, options);

    y(:,1) = y0;
    j = 2;
    
    for i = 2:length(t)
        y0 = y0 + (t(i) - t(i-1)) * odefun(t(i-1), y0);

        if i == saveIndices(j)
            y(:,j) = y0;
            j = j + 1;
        end
    end

    [tOut, y] = functionOutputParser(tOut, y, nargout);
end
