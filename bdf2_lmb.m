function y = bdf2_lmb(odefun,t,y0)
    verbose = false;
    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;
    y(:,2) = lmgb(@(h) h - (y0 + (t(2)-t(1))*odefun(t(2),h)),y0);
    
    for i = 3:n
        [y(:,i),fval,exitflag] = lmgb(@(h) h - (4/3 * y(:,i-1) - 1/3 * y(:,i-2) + 2/3 * (t(i)-t(i-1)) * odefun(t(i), h)),y(:,i-1));
        if verbose == true
            fprintf('fval: %i, exitflag: %g\n',norm(fval,'inf'),exitflag);
        end
    end
end