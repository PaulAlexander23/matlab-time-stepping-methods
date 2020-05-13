function [tOut, y] = bdf2(odefun,tOut,yn,options)
    if nargin < 4
        options = odeset();
    end
    options = ensureSolverSet(options);

    n = length(tOut);
    y = zeros(length(yn),n);
    
    [t, saveIndices] = timepointsWithMaxStep(tOut, options);
    validateTimeStepsEqual(t);

    y(:,1) = yn;
    j = 2;

    yCoeff1 = [1, -1]';
    fCoeff1 = -1;
    yCoeff2 = [1, -4/3, 1/3]';
    fCoeff2 = -2/3;

    for i = 2:length(t)
        dt = t(i)-t(i-1);
        if i == 2
            explicitF = yn * yCoeff1(2);

            ynm1 = yn;
            yn = options.optimmethod( ...
                @(h) funBDF1(odefun, t(i), h, dt, explicitF, options), ...
                yn, ...
                options.optimoptions);
        else
            explicitF = [yn, ynm1] * yCoeff2(2:end);

            ynm1 = yn;
            yn = options.optimmethod( ...
                @(h) funBDF2(odefun, t(i), h, dt, explicitF, options), ...
                yn, ...
                options.optimoptions);
        end

        if i == saveIndices(j)
            y(:,j) = yn;
            j = j + 1;
        end

        if any(isnan(yn))
            fprintf('Nan`s in solution\n')
            break;
        end
    end

    function [F,J] = funBDF1(odefun, t, h, dt, explicitF, options)
        f = odefun(t,h);
        F = h * yCoeff1(1) + dt * f * fCoeff1(1) + explicitF;

        if nargout == 2
            jac = options.Jacobian(t,h);
            J = speye(length(h)) * yCoeff1(1) + dt * jac * fCoeff1(1);
        end
    end

    function [F,J] = funBDF2(odefun, t, h, dt, explicitF, options)
        f = odefun(t,h);
        F = h * yCoeff2(1) + dt * f * fCoeff2(1) + explicitF;

        if nargout == 2
            jac = options.Jacobian(t,h);
            J = speye(length(h)) * yCoeff2(1) + dt * jac * fCoeff2(1);
        end
    end
    
    [tOut, y] = functionOutputParser(tOut, y, nargout);
end
