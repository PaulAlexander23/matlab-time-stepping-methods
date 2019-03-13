function y = am3(odefun,t,y0)
    n = length(t);
    y = zeros(length(y0), n);
    
    y(:,1) = y0;
    y(:,2) = y(:,1);
    m = 1;
    err = 1;
    while err > 1e-6 && m <= 20
        temp = y(:,2);
        y(:,2) = y(:,1) + (t(2)-t(1))/2 * (odefun(t(2), y(:,2)) + odefun(t(1), y(:,1)));
        err = norm(y(:,2) - temp)/norm(temp);
        m = m + 1;
    end
    
    for i = 3:n
        y(:,i) = y(:,i-1);
        dt = t(i) - t(i-1);
        m = 1;
        err = 1;
        while err > 1e-6 && m <= 20
            temp = y(:,i);
            y(:,i) = y(:,i-1) + dt/12 * (5 * odefun(t(i), y(:,i)) + 8 * odefun(t(i-1), y(:,i-1)) - odefun(t(i-2), y(:,i-2)));
            err = norm(y(:,i) - temp)/norm(temp);
            m = m + 1;
        end
    end
end