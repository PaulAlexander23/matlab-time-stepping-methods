function [t, y] = abm1(odefun,t,y0,options)
    if nargin < 4, options = struct(); end

    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;
    
    for i = 2:n
        dt = t(i)-t(i-1);
        y(:,i) =  y(:,i-1) + dt * odefun(t(i-1), y(:,i-1)); % AB1
        
        m = 1;
        err = 1;
        while err > 1e-6 && m <= 20
            temp = y(:,i);
            y(:,i) = y(:,i-1) + dt * odefun(t(i), y(:,i)); % AM1
            err = norm(y(:,i) - temp)/norm(temp);
            m = m + 1;
        end
    end

    [t, y] = functionOutputParser(t, y, nargout);
end
