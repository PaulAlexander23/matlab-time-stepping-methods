function y = bdf3_semi_implicit(odefun, t, y0, optimmethod, explicitfun)
    n = length(t);
    y = zeros(length(y0), n);

    y(:, 1:3) = bdf2_semi_implicit(odefun, t(1:3), y0, optimmethod, explicitfun);

    function [F, J] = fun(h, odefun, i, t, y, expf)

        if nargout == 1
            f = odefun(t, h);
            F = (h - (18 * y(:, i - 1) - 9 * y(:, i - 2) + ...
                2 * y(:, i - 3)) / 11) / (t(i) - t(i - 1)) - 6/11 * f;
        elseif nargout == 2
            [f, j] = odefun(t, h);
            F = (h - (18 * y(:, i - 1) - 9 * y(:, i - 2) + ...
                2 * y(:, i - 3)) / 11) / (t(i) - t(i - 1)) - 6/11 * f;
            J = speye(length(f)) / (t(i) - t(i - 1)) - 2/3 * j;
        end

        F = F - expf;
        
    end

    for i = 4:n
        expf = 3 * explicitfun(t, y(:, i - 1)) - ...
            3 * explicitfun(t, y(:, i - 2)) + ...
            explicitfun(t, y(:, i - 3));
        y(:, i) = optimmethod(@(h) fun(h, odefun, i, t, y, expf), ...
            y(:, i - 1));

        if any(isnan(y(:, i)))
            fprintf('Nan`s in solution\n')
            break;
        end

    end

end
