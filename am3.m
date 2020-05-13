function [tOut, y] = am3(odefun,tOut,yn,options)
    if nargin < 4, options = odeset(); end

    n = length(tOut);
    y = zeros(length(yn), n);
    
    [t, saveIndices] = timepointsWithMaxStep(tOut, options);
    validateTimeStepsEqual(t);
    
    y(:,1) = yn;
    j = 2;

    for i = 2:length(t)
        dt = t(i) - t(i-1);
        m = 1;
        err = 1;
        
        if i == 2
            ynm1 = yn;
            while err > 1e-6 && m <= 20
                temp = yn;
                yn = ynm1 + (t(2)-t(1))/2 * (odefun(t(2), yn) + odefun(t(1), ynm1));
                err = norm(yn - temp)/norm(temp);
                m = m + 1;
            end
        else
            ynm2 = ynm1;
            ynm1 = yn;
            while err > 1e-6 && m <= 20
                temp = yn;
                yn = ynm1 + dt/12 * (5 * odefun(t(i), yn) + 8 * odefun(t(i-1), ynm1) - odefun(t(i-2), ynm2));
                err = norm(yn - temp)/norm(temp);
                m = m + 1;
            end
        end

        if i == saveIndices(j)
            y(:,j) = yn;
            j = j + 1;
        end
    end

    [tOut, y] = functionOutputParser(tOut, y, nargout);
end
