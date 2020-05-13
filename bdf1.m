function [tOut, y] = bdf1(odefun,tOut,yn,options)
    if nargin < 4
        options = odeset();
    end
    options = ensureSolverSet(options);

    n = length(tOut);
    y = zeros(length(yn),n);
    
    [t, saveIndices] = timepointsWithMaxStep(tOut, options);

    y(:,1) = yn;
    j = 2;
    
    yCoeff = [1, -1]';
    fCoeff = -1;

    for i = 2:length(t)
        dt = t(i) - t(i-1);
        explicitF = yn * yCoeff(2);

        yn = options.optimmethod( ...
            @(h) fun(odefun, t(i), h, dt, explicitF, options), ...
            yn, ...
            options.optimoptions);
        
        if i == saveIndices(j)
            y(:,j) = yn;
            j = j + 1;
        end

        if any(isnan(yn))
            fprintf('Nan`s in solution\n')
            break;
        end
    end

    function [F,J] = fun(odefun, t, h, dt, explicitF, options)
        f = odefun(t,h);
        F = h * yCoeff(1) + dt * f * fCoeff(1) + explicitF;

        if nargout == 2
            jac = options.Jacobian(t,h);
            J = speye(length(h)) * yCoeff(1) + dt * jac * fCoeff(1);
        end
    end
    
    [tOut, y] = functionOutputParser(tOut, y, nargout);
end
