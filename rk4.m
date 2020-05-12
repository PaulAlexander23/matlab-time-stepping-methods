function [tOut, y] = rk4(odefun, tOut, yn, options)
    if nargin < 4, options = odeset(); end

    n = length(tOut);
    y = zeros(length(yn),n);
    
    [t, saveIndices] = timepointsWithMaxStep(tOut, options);

    y(:,1) = yn;
    j = 2;
    
    for i = 2:length(t)
        h = t(i)-t(i-1);
        k1 = h * odefun(t(i-1), yn);
        k2 = h * odefun(t(i-1) + h/2, yn + k1/2);
        k3 = h * odefun(t(i-1) + h/2, yn + k2/2);
        k4 = h * odefun(t(i-1) + h, yn + k3);
        
        yn = yn + (k1 + 2 * k2 + 2 * k3 + k4)/6;

        if i == saveIndices(j)
            y(:,j) = yn;
            j = j + 1;
        end
    end
    
    [tOut, y] = functionOutputParser(tOut, y, nargout);
end
