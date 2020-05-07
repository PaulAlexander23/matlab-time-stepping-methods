function [t, y] = bdf1(odefun,t,y0,options)
    if nargin < 4
        options = odeset();
    end
    options = ensureSolverSet(options);

    n = length(t);
    y = zeros(length(y0),n);
    
    y(:,1) = y0;
    
    yCoeff = [1, -1]';
    fCoeff = -1;

    for i = 2:n
        dt = t(i) - t(i-1);

        explicitF = y(:,i-1) * yCoeff(2);

        y(:,i) = options.optimmethod( ...
            @(h) fun(odefun, t(i), h, dt, explicitF, options), ...
            y(:,i-1), ...
            options.optimoptions);
        
        if any(isnan(y(:,i)))
            fprintf('Nan`s in solution\n')
            break;
        end
    end

    function [F,J] = fun(odefun, t, h, dt, explicitF, options)
        f = odefun(t,h);
        F = h * yCoeff(1) + dt * f * fCoeff(1) + explicitF;

        if nargout == 2
            j = options.Jacobian(t,h);
            J = speye(length(h)) * yCoeff(1) + dt * j * fCoeff(1);
        end
    end
    
    [t, y] = functionOutputParser(t, y, nargout);
end
