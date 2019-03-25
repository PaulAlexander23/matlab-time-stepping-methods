function y = bdf1_lmb_flux(fluxfun,t,y0,x)
    verbose = false;
    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;

    for i = 2:n
        [y(:,i),fval,exitflag] = lmgb(@(h) h - (y(:,i-1) + (t(i)-t(i-1)) * div(x,fluxfun(t, h))),y(:,i-1));
        if verbose == true
            fprintf('fval: %i, exitflag: %g\n',norm(fval,'inf'),exitflag);
        end
        if any(isnan(y(:,i)),'all')
            fprintf('Nan`s in solution\n')
            break;
        end
    end
    
    function out = div(x,y)
        divdeg = cat(3,[1,0;0,0],[0,0;0,1]);
        d2y = diff_ps_2d(x,y,divdeg);
        out = d2y(:,:,1) + d2y(:,:,2);
    end
end