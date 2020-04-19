function [t, y] = ab1(odefun,t,y0)
    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;
    
    for i = 2:n
        y(:,i) = y(:,i-1) + (t(i)-t(i-1)) * odefun(t(i-1), y(:,i-1));
    end

    [t, y] = functionOutputParser(t, y, nargout);
end
