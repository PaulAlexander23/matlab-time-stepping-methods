function [y,te,ye,ie] = bdf1(odefun,t,y0,optimmethod,event)
    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;
    if nargin > 4
        eValue = event(0,y0);
        ie = [];
        te = [];
        ye = [];
    end
    
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
        y(:,i) = optimmethod(@(h) fun(h,odefun,i,t,y),y(:,i-1));
        
        if any(isnan(y(:,i)))
            fprintf('Nan`s in solution\n')
            break;
        end
        if nargin > 4
            temp = eValue;
            [eValue,eTerminal,eDirection] = event(t(i),y(:,i));
            eIndex = (temp .* eValue < 0) .* (eDirection .* temp <= 0);
            if any(ie)
                ie = [ie,eIndex];
                te = [te,t(i)];
                ye = [ye,y(:,i)];
            end
            if any(eTerminal(eIndex),'all')
                break;
            end
        end
    end
end