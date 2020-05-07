function [t, y] = bdf2(odefun,t,y0,options)
    if nargin < 4
        options = odeset();
    end
    options = ensureSolverSet(options);

    n = length(t);
    y = zeros(length(y0),n);
    
    [~, y0] = bdf1(odefun,t(1:2),y0,options);
    y(:,1:2) = y0';

    function [F,J] = fun(h,odefun,i,t,y)
        if nargout == 1
            f = odefun(t,h);
            F = (h - (4*y(:,i-1)-y(:,i-2))/3)/(t(i)-t(i-1)) - 2/3*f;
        elseif nargout == 2
            [f, j] = odefun(t,h);
            F = (h - (4*y(:,i-1)-y(:,i-2))/3)/(t(i)-t(i-1)) - 2/3*f;
            J = speye(length(f))/(t(i)-t(i-1)) - 2/3*j;
        end
    end
    
    for i = 3:n
        y(:,i) = options.optimmethod( ...
            @(h) fun(h,odefun,i,t,y), ...
            y(:,i-1),...
            options.optimoptions);
        if any(isnan(y(:,i)))
            fprintf('Nan`s in solution\n')
            break;
        end
    end

    [t, y] = functionOutputParser(t, y, nargout);
end
