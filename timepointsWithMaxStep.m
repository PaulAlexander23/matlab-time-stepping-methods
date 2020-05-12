function [t, saveIndices] = timepointsWithMaxStep(tOut, options)
    t = tOut;
    if isfield(options, 'MaxStep')
        if ~isempty(options.MaxStep)
            t = cell2mat(arrayfun(@(x,y)(x:options.MaxStep:y)',...
                tOut, [tOut(2:end); tOut(end)], ...
                'UniformOutput',false));
        end
    end

    t = unique(t);
    
    [saveIndices, ~] = find(tOut' == t);
end
