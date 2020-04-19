function [t, y] = bdf3(odefun,t,y0,options)
    if nargin < 4
        options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
            optimoptions('fsolve', 'Display', 'off')));
    end
    n = length(t);
    y = zeros(length(y0),n);
    
    [~, y(:,1:3)] = bdf2(odefun,t(1:3),y0,options);
 
    function [F,J] = fun(h,odefun,i,t,y)
        if nargout == 1
            f = odefun(t,h);
            F = (h - (18*y(:,i-1)-9*y(:,i-2)+2*y(:,i-3))/11)/(t(i)-t(i-1)) - 6/11*f;
        elseif nargout == 2
            [f, j] = odefun(t,h);
            F = (h - (18*y(:,i-1)-9*y(:,i-2)+2*y(:,i-3))/11)/(t(i)-t(i-1)) - 6/11*f;
            J = speye(length(f))/(t(i)-t(i-1)) - 2/3*j;
        end
    end
    
    for i = 4:n
        y(:,i) = options.optimmethod(@(h) fun(h,odefun,i,t,y),y(:,i-1));
        if any(isnan(y(:,i)))
            fprintf('Nan`s in solution\n')
            break;
        end
    end

    [t, y] = functionOutputParser(t, y, nargout);
end
