function y = bdf1_newton_jacobian(odefun,t,y0)
    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;
    
    function [F,J] = fun(h,odefun,i,t,y)
        [f, j] = odefun(t,h);
        F = (h - y(:,i-1))/(t(i)-t(i-1)) - f;
        J = 1/(t(i)-t(i-1)) - j;
    end
    
    for i = 2:n
        y(:,i) = newton(@(h) fun(h,odefun,i,t,y),y(:,i-1));
        if any(isnan(y(:,i)),'all')
            fprintf('Nan`s in solution\n')
            break;
        end
    end
end