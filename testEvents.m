function tests = testEvents()
    tests = functiontests(localfunctions);
end

function testEventODE45(testCase)
    [t, y0, options] = setupVars();

    sol = ode45(@odefun, t, y0, options);

    verifyZeroAfterEvent(testCase, sol);
end

function testEventBDF1SI(testCase)
    [t, y0, options] = setupVars();

    sol = bdf1si(odefunSemiImplicit(), t, y0, options);

    verifyZeroAfterEvent(testCase, sol);
end

function testEventBDF2SI(testCase)
    [t, y0, options] = setupVars();

    sol = bdf2si(odefunSemiImplicit(), t, y0, options);

    verifyZeroAfterEvent(testCase, sol);
end

function testEventBDF3SI(testCase)
    [t, y0, options] = setupVars();

    sol = bdf3si(odefunSemiImplicit(), t, y0, options);

    verifyZeroAfterEvent(testCase, sol);
end

function testEventBDF4SI(testCase)
    [t, y0, options] = setupVars();

    sol = bdf4si(odefunSemiImplicit(), t, y0, options);

    verifyZeroAfterEvent(testCase, sol);
end

function testEventBDF5SI(testCase)
    [t, y0, options] = setupVars();

    sol = bdf5si(odefunSemiImplicit(), t, y0, options);

    verifyZeroAfterEvent(testCase, sol);
end

function testEventBDF6SI(testCase)
    [t, y0, options] = setupVars();

    sol = bdf6si(odefunSemiImplicit(), t, y0, options);

    verifyZeroAfterEvent(testCase, sol);
end

% End of test functions

function [t, y0, options] = setupVars()
    t = [0,0.5,1.5,2]';
    y0 = 1;
    options = odeset("event",@myEvent, 'MaxStep', 1e-2);
end

function f = odefun(t, y)
    f = -1 + 0*y + 0*t;
end

function f = odefunSemiImplicit()
    f = struct();
    f.explicit = @(t, y) -1/2 + 0*y + 0*t;
    f.implicit = @(t, y) -1/2 + 0*y + 0*t;
end

function [value, isterminal, direction] = myEvent(~, y)
    value = y;
    isterminal = true;
    direction = 0;
end

function verifyZeroAfterEvent(testCase, sol)
    index = sol.x > 1 + 1e-6;
    verifyEqual(testCase, sol.y(index), 0*sol.y(index), 'AbsTol', 1e-6)
end
