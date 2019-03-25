function y = bdf1_newton(odefun,t,y0)
    verbose = true;
    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;
    if verbose == false
        options = optimoptions('fsolve','Display','none');
    else
        options = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-5*numel(y0));
    end
    
    for i = 2:n
        y(:,i) = fsolve(@(h) h - (y(:,i-1) + (t(i)-t(i-1)) * odefun(t, h)),y(:,i-1),options);
        if any(isnan(y(:,i)),'all')
            fprintf('Nan`s in solution\n')
            break;
        end
    end
end