function y = abm2(odefun,t,y0)
    n = length(t);
    y = zeros([size(y0),n]);
    
    y(:,:,1) = y0;
    
    dt = t(2)-t(1);
    y(:,:,2) = y(:,:,1) + dt * odefun(t(1), y(:,:,1)); % AB1
    m = 1;
    err = 1;
    while err > 1e-6 && m <= 20
        temp = y(:,:,2);
        y(:,:,2) = y(:,:,1) + dt/2 * (odefun(t(2), y(:,:,2)) + odefun(t(1), y(:,:,1))); % AM1
        err = norm(y(:,:,2) - temp)/norm(temp);
        m = m + 1;
    end
    
    for i = 3:n
        dt = t(i)-t(i-1);
        y(:,:,i) =  y(:,:,i-1) + dt/2 * (3 * odefun(t(i-1), y(:,:,i-1)) - odefun(t(i-2), y(:,:,2))); % AB2
        m = 1;
        err = 1;
        while err > 1e-6 && m <= 20
            temp = y(:,:,i);
            y(:,:,i) = y(:,:,i-1) + dt/2 * (odefun(t(i), y(:,:,i)) + odefun(t(i-1), y(:,:,i-1))); % AM2
            err = norm(y(:,:,i) - temp)/norm(temp);
            m = m + 1;
        end
    end
    
    y = squeeze(y);
    
end