function [quit, value, ie, xe, ye] = handleEvents(event, t, y, oldValue)

    quit = false;
    ie = 0;
    xe = 0;
    ye = 0*y;

    [value, isterminal, direction] = event(t, y);

    signChanged = (sign(value) - sign(oldValue))/2;

    if any(signChanged .* ~direction) || ...
            any(abs(signChanged) .* (signChanged == direction))
        ie = find(signChanged);
        xe = t;
        ye = y;

        if any(isterminal(ie))
            quit = true;
        end
    end
end
