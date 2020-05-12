function validateTimeStepsEqual(t)
    if any(abs(diff(diff(t))) > eps)
        error(['Step sizes must be equal. First 10 steps:', ...
            sprintf('%g, ',diff(t(1:11)))]);
    end
end
