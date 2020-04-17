function y = ab3cn(odefunNL, t, y0, options)
    n = length(t);
    y = zeros(length(y0),n);
    
    y(:, 1) = y0;
    y(:, 2) = (1 - (t(2) - t(1)) * options.Linear(t(2))) \ ...
        (y(:, 1) + (t(2) - t(1)) * odefunNL(t(1), y(:, 1)));
    y(:, 3) = (1 - (t(3) - t(2)) * options.Linear(t(3))) \ ...
        (y(:, 2) + (t(3) - t(2)) * (3/2 * odefunNL(t(2), y(:, 2)) - 1/2 * odefunNL(t(1), y(:, 1))));
    
    coeff = [23/12,-4/3,5/12]';
    
    for i = 4:n
        tau = t(i) - t(i - 1);
        y(:, i) = (1 - tau/2 * options.Linear(t(i))) \ ...
            (y(:, i-1) + ...
            tau * odefunNL(t(i-1:-1:i-3),y(:,i-1:-1:i-3)) * coeff + ...
            tau/2 * options.Linear(t(i)) * y(:, i-1));
    end
end