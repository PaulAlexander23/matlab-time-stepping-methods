function y = bdf3_lmb(odefun,t,y0)
    verbose = false;
    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;
    y(:,2) = lmgb(@(h) h - (y0 + (t(2)-t(1))*odefun(t(2),h)),y(:,1));
    y(:,3) = lmgb(@(h) h - (4/3 * y(:,2) - 1/3 * y(:,1) + 2/3 * (t(3)-t(2)) * odefun(t(3), h)),y(:,2));
    
    for i = 4:n
        [y(:,i),fval,exitflag] = lmgb(@(h) h - (18 * y(:,i-1) - 9 * y(:,i-2) + 2 * y(:,i-3) + 6 * (t(i)-t(i-1)) * odefun(t(i), h))/11,y(:,i-1));
        if verbose == true
            fprintf('fval: %i, exitflag: %g\n',norm(fval,'inf'),exitflag);
        end
    end
end