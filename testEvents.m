function tests = testEvents()
    tests = functiontests(localfunctions);
end

function testEvent(testCase)
    odefun = @(t, y) -1;
    t = linspace(0,2,3)';
    y0 = 1;
    options = odeset("event",@myEvent, 'MaxStep', 1e-3);

    sol = ode45(odefun, t, y0, options);

    save("temp.mat")

    function [value, isterminal, direction] = myEvent(t, y)
        value = y;
        isterminal = true;
        direction = 0;
    end
end

function testEventBDF2SI(testCase)
    explicitOdefun = @(t,y) -1/2;
    implicitOdefun = @(t,y) -1/2;
    odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
    t = [0,0.5,1.5,2]';
    y0 = 1;
    options = odeset("event",@myEvent, 'MaxStep', 1e-3);

    sol = bdf2si(odefun, t, y0, options);

    save("temp.mat")

    function [value, isterminal, direction] = myEvent(t, y)
        value = y;
        isterminal = true;
        direction = 0;
    end
end
