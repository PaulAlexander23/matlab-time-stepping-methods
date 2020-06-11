function [t, y] = functionOutputParserEvents(t, y, ie, xe, ye, passedNargout)
    if passedNargout == 1
        t = struct('x', t, 'y', y.', 'ie', ie, 'xe', xe, 'ye', ye.');
    elseif passedNargout == 2
        y = y.';
    end
end

