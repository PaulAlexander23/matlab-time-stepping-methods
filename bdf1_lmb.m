function y = bdf1_lmb(odefun,t,y0)
    verbose = false;
    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;

    for i = 2:n
        [y(:,i),fval,exitflag] = lmgb(@(h) h - (y(:,i-1) + (t(i)-t(i-1)) * odefun(t, h)),y(:,i-1));
        if verbose == true
            fprintf('fval: %i, exitflag: %g\n',norm(fval,'inf'),exitflag);
        end
    end
end