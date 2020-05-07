function options = ensureSolverSet(options)
    if ~isfield(options, 'optimmethod')
        options.optimmethod = @fsolve;
    end
    if ~isfield(options, 'optimoptions')
        options.optimoptions = optimoptions('fsolve', 'Display', 'off');
    end

end
