function [t, y] = functionOutputParser(t, y, passedNargout)
    if passedNargout == 1
        t = struct('x', t, 'y', y');
    elseif passedNargout == 2
        y = y';
    end
end
