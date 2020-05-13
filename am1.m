function [tOut, y] = am1(odefun,tOut,yn,options)
    if nargin < 4, options = odeset(); end

    n = length(tOut);
    y = zeros(length(yn), n);
    
    [t, saveIndices] = timepointsWithMaxStep(tOut, options);

    y(:,1) = yn;
    j = 2;

    for i = 2:length(t)
        ynm1 = yn;
        dt = t(i) - t(i-1);

        m = 1;
        err = 1;
        while err > 1e-6 && m <= 20
            temp = yn;
            yn = ynm1 + dt * odefun(t(i), yn);
            err = norm(yn - temp)/norm(temp);
            m = m + 1;
        end

        if i == saveIndices(j)
            y(:,j) = yn;
            j = j + 1;
        end
    end

    [tOut, y] = functionOutputParser(tOut, y, nargout);
end
