function y = bdf(odefun, t, y0, options)
    %BDF Solve differential equations using backward differentiation formulae.
    % Handle options
    if nargin < 4, options = struct(); end
    
    if ~isfield(options,'Order'), options.Order = 4; end
    if ~isfield(options,'Optimmethod')
        optimopts = optimoptions('fsolve','Display','off');
        options.Optimmethod = @(fun, x0) fsolve(fun, x0, optimopts);
    end
    if ~isfield(options,'Jacobian'), options.Jacobian = []; end
    
    % Preallocate solution
    n = length(t);
    y = zeros(length(y0), n);
    
    % Initialise first steps iteratively with lower orders
    if options.Order == 1
        y(:, 1) = y0;
    else
        newoptions = options;
        newoptions.Order = options.Order - 1;
        y(:, 1:options.Order) = bdf(odefun, t(1:options.Order), ...
            y0, newoptions);
    end
    
    % Compute BDF coefficients
    m = (1:options.Order);
    k = triu(repmat((0:options.Order)',1,options.Order),-1);
    bin = factorial(m)./(factorial(m-k).*factorial(k));
    bin = triu(bin,-1).* (-1).^(0:options.Order)';
    v = 1./m;
    coeff = bin*v';
    
    % Compute interpolation coefficients
    interpCoeff = -bin(2:end,end);
    
    % Define function to solve
    function [F,J] = fun(h,odefun,t,dt,coeff,expf,options)
        % BDF
        F = coeff * h - odefun(t, h) * dt - expf;
        
        if ~isempty(options.Jacobian)
            % Generate Jacobian
            j = options.Jacobian(t, h);
            % Compute BDF jacobian
            J = coeff * speye(size(j)) - j * dt;
        end
    end
    
    % Time stepping loop
    for i = options.Order+1:n
        % BDF constant part
        expf = - y(:, i-1:-1:i-options.Order) * coeff(2:end);
        
        % Semi implicit function interpolation
        if isfield(options,'ExplicitFcn')
            expfeval = options.ExplicitFcn(t(i-1:-1:i-options.Order), ...
                y(:, i-1:-1:i-options.Order));
            expf = expf + (t(i) - t(i-1)) * expfeval * interpCoeff;
        end
        % BDF implicit solve
        y(:,i) = options.Optimmethod( ...
            @(h) fun(h, odefun, t(i), t(i) - t(i - 1), coeff(1), expf, ...
            options), y(:, i - 1));
        
        % Output handling
        if any(isnan(y(:,i)))
            fprintf('Nan`s in solution\n')
            break;
        end
    end
end