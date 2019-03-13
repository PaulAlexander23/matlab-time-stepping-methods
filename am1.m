function y = am1(odefun,t,y0)
    n = length(t);
    y = zeros(length(y0), n);
    
    y(:,1) = y0;
    
    for i = 2:n
        y(:,i) = y(:,i-1);
        dt = t(i) - t(i-1);
        m = 1;
        err = 1;
        while err > 1e-6 && m <= 20
            temp = y(:,i);
            y(:,i) = y(:,i-1) + dt * odefun(t(i), y(:,i));
            err = norm(y(:,i) - temp)/norm(temp);
            m = m + 1;
        end
    end
end