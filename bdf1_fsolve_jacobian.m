function y = bdf1_newton_jacobian(odefun,t,y0,jacobian)
    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;
    options = optimoptions('fsolve','SpecifyObjectiveGradient',true);

    function [F,J] = fun(h,odefun,i,t,y,jacobian)
        F = h - (y(:,i-1) + (t(i)-t(i-1)) * odefun(t, h));
        if nargout > 1
            J = 1 - (t(i)-t(i-1))*jacobian(h);
        end
    end
    for i = 2:n
        y(:,i) = fsolve(@(h) fun(h,odefun,i,t,y,jacobian),y(:,i-1),options);
        if any(isnan(y(:,i)),'all')
            fprintf('Nan`s in solution\n')
            break;
        end
    end
end