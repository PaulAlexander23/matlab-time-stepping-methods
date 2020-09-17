function validateTimeStepsEqual(t)
    if any(abs(diff(diff(t))) > 1e-12)
        error(['Step sizes must be equal. ', ...
            'Largest offset: ', sprintf('%g. ', max(abs(diff(diff(t))))), ...
            'First 10 step sizes: ', sprintf('%g, ',diff(t(1:11)))]);
    end
end
