function y = bdf1_newton_flux(fluxfun,t,y0,x,fluxjacobian,method)
    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;
    options = optimoptions('fsolve','SpecifyObjectiveGradient',true);
    
    for i = 2:n
        y(:,i) = fsolve(@(t,h) fun(h,fluxfun,i,t,y,x,fluxjacobian,method),...
            y(:,i-1),options);
        if any(isnan(y(:,i)),'all')
            fprintf('Nan`s in solution\n')
            break;
        end
    end
    
    function [f,J] = fun(h,fluxfun,i,t,y,x,fluxjacobian,method)
        f = (h - y(:,i-1))/ (t(i)-t(i-1)) + div(x,fluxfun(t, h),method);
        if nargout > 1
            J = 1 / (t(i)-t(i-1)) - div(x,fluxjacobian(h),method);
        end
    end
    
    function out = div(x,y,method)
        out = cell2mat(method(x,y(:,:,1),[1,0]')) + ...
            cell2mat(method(x,y(:,:,2),[0,1]'));
    end
end