function y = bdf1_semi_implicit(odefun, t, y0, optimmethod, explicitfun)
    n = length(t);
    y = zeros(length(y0), n);

    y(:, 1) = y0;

    function [F, J] = fun(h, odefun, i, t, y, expf)

        if nargout == 1
            f = odefun(t(i), h);
            F = (h - y(:, i - 1)) / (t(i) - t(i - 1)) - f;
        elseif nargout == 2
            [f, j] = odefun(t(i), h);
            F = (h - y(:, i - 1)) / (t(i) - t(i - 1)) - f;
            J = speye(length(f)) / (t(i) - t(i - 1)) - j;
        end
        
        F = F - expf;
        
    end

    for i = 2:n
        expf = explicitfun(t(i), y(:, i - 1));
        y(:, i) = optimmethod(@(h) fun(h, odefun, i, t, y, expf), ...
            y(:, i - 1));

        if any(isnan(y(:, i)))
            fprintf('Nan`s in solution\n')
            break;
        end

    end

end
