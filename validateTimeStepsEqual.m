function validateTimeStepsEqual(t)
    if any(abs(diff(diff(t))) > 1e-13)
        error(['Step sizes must be equal. First 10 step sizes:', ...
            sprintf('%g, ',diff(t(1:11))), '\n']);
    end
end
