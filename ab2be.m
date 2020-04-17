function y = ab2be(odefunLinear, odefunNonlinear, t, y0, optimmethod)
    n = length(t);
    y = zeros(length(y0),n);
    
    %Change wind up
    odefun = @(t, y) odefunLinear(t, y) + odefunNonlinear(t, y);
    y(:,1) = y0;
    y(:,2) = optimmethod(@(x) x - y(:, 1) - ...
        (t(2) - t(1)) * odefunNonlinear(t(1), y(:, 1)) - ...
        (t(2) - t(1)) * odefunLinear(t(2), x), ...
        y(:, 1));
    
    coeff = [3/2,-1/2]';
    
    for i = 3:n
        tau = t(i) - t(i - 1);
        y(:, i) = optimmethod(@(x) x - y(:, i-1) - ...
            tau * odefunNonlinear(t(i-1:-1:i-2), y(:, i-1:-1:i-2)) * coeff - ...
            tau * odefunLinear(t(i), x), ...
            y(:, i - 1));
    end
end