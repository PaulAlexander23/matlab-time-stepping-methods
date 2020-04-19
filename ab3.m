function [t, y] = ab3(odefun,t,y0)
    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;
    y(:,2) = y(:,1) + (t(2)-t(1)) * odefun(t(1), y(:,1));
    y(:,3) = y(:,2) + (t(3)-t(2)) * (3/2 * odefun(t(2), y(:,2)) - 1/2 * odefun(t(1), y(:,1)));
    
    coeff = [23/12,-4/3,5/12]';
    
    for i = 4:n
        y(:,i) = y(:,i-1) + (t(i)-t(i-1)) * odefun(t(i-1:-1:i-3),y(:,i-1:-1:i-3)) * coeff;
    end

    [t, y] = functionOutputParser(t, y, nargout);
end
