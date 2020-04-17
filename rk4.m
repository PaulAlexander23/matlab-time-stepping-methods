function y = rk4(odefun,t,y0)
    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;
    
    for i = 1:n-1
        h = t(i+1)-t(i);
        
        k1 = h * odefun(t(i), y(:,i));
        k2 = h * odefun(t(i) + h/2, y(:,i) + k1/2);
        k3 = h * odefun(t(i) + h/2, y(:,i) + k2/2);
        k4 = h * odefun(t(i) + h, y(:,i) + k3);
        
        y(:,i+1) = y(:,i) + (k1 + 2 * k2 + 2 * k3 + k4)/6;
    end
    
    y = y.';
end
