function [t, y] = abm3(odefun,t,y0,options)
    if nargin < 4, options = odeset(); end

    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;
    
    dt = t(2)-t(1);
    y(:,2) = y(:,1) + dt * odefun(t(1), y(:,1)); % AB1
    m = 1;
    err = 1;
    while err > 1e-6 && m <= 20
        temp = y(:,2);
        y(:,2) = y(:,1) + dt/2 * (odefun(t(2), y(:,2)) + odefun(t(1), y(:,1))); % AM1
        err = norm(y(:,2) - temp)/norm(temp);
        m = m + 1;
    end
    
    dt = t(3)-t(2);
    y(:,3) =  y(:,2) + dt/2 * (3 * odefun(t(2), y(:,2)) - odefun(t(1), y(:,1))); % AB2
    m = 1;
    err = 1;
    while err > 1e-6 && m <= 20
        temp = y(:,3);
        y(:,3) = y(:,2) + dt/2 * (odefun(t(3), y(:,3)) + odefun(t(2), y(:,2))); % AM2
        err = norm(y(:,3) - temp)/norm(temp);
        m = m + 1;
    end
    
    coeff = [23/12,-4/3,5/12]';
    
    for i = 4:n
        dt = t(i)-t(i-1);
        y(:,i) = y(:,i-1) + dt * odefun(t(i-1:-1:i-3),y(:,i-1:-1:i-3)) * coeff; % AB3
        m = 1;
        err = 1;
        while err > 1e-6 && m <= 20
            temp = y(:,i);
            y(:,i) = y(:,i-1) + dt/12 * (5 * odefun(t(i), y(:,i)) + 8 * odefun(t(i-1), y(:,i-1)) - odefun(t(i-2), y(:,i-2))); % AM3
            err = norm(y(:,i) - temp)/norm(temp);
            m = m + 1;
        end
    end
    
    [t, y] = functionOutputParser(t, y, nargout);
end
