function y = ab2(odefun,t,y0)
    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;
    y(:,2) = y0 + (t(2)-t(1)) * odefun(t(1), y0);
    
    for i = 3:n
        y(:,i) = y(:,i-1) + (t(i)-t(i-1)) * (3/2 * odefun(t(i-1), y(:,i-1)) - 1/2 * odefun(t(i-2), y(:,i-2)));
    end
end