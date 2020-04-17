function y = bdf1(odefun,t,y0,options)
    if nargin < 4
        options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
            optimoptions('fsolve', 'Display', 'off')));
    end

    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;
    
    function [F,J] = fun(h,odefun,i,t,y)
        if nargout == 1
            f = odefun(t,h);
            F = (h - y(:,i-1))/(t(i)-t(i-1)) - f;
        elseif nargout == 2
            [f, j] = odefun(t,h);
            F = (h - y(:,i-1))/(t(i)-t(i-1)) - f;
            J = speye(length(f))/(t(i)-t(i-1)) - j;
        end
    end
    
    for i = 2:n
        y(:,i) = options.optimmethod(@(h) fun(h,odefun,i,t,y),y(:,i-1));
        
        if any(isnan(y(:,i)))
            fprintf('Nan`s in solution\n')
            break;
        end
    end
end
