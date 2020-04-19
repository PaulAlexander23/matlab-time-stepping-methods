function tests = testSemiImplicitSolvers()
    tests = functiontests(localfunctions);
end

function testOdefunSolve(testCase)
    explicitOdefun = @(t,y) y;
    implicitOdefun = @(t,y) -y^2;
    t = linspace(0,10,100)';
    y0 = 0.1;
    options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
        optimoptions('fsolve', 'Display', 'off')));

    solverList = {@ab1be, @ab2be, @ab3cn, @bdf1si, @bdf2si, @bdf3si};
    expectedAccuracy = {1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8};

    A = y0 / (1 - y0);
    expected = A*exp(t')./(1 + A*exp(t'));
    plot(expected); hold on;
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        solution = solverList{n}(explicitOdefun, implicitOdefun, t, y0, options);
        actual = solution;
        plot(actual);

        verifyEqual(testCase, actual, expected, 'AbsTol', expectedAccuracy{n});
    end
end

function testConvergenceRates(testCase)
    explicitOdefun = @(t,y) y;
    implicitOdefun = @(t,y) -y^2;
    tN = round(logspace(2,2.5,6)); %[100, 110, 120];
    tL = 1;
    y0 = 0.5;
    options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
        optimoptions('fsolve', 'Display', 'off')));

    solverList = {@ab1be, @ab2be, @ab3cn, @bdf1si, @bdf2si, @bdf3si};
    expected = {2, 2, 3, 2, 3, 4};

    difference = zeros(length(tN),1);
    trueValue = exp(tL)./(1 + exp(tL));
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));

        for m = 1:length(tN)
            t = linspace(0,tL,tN(m))';
            solution = solverList{n}(explicitOdefun, implicitOdefun, t, y0, options);
            difference(m) = solution(end) - trueValue;
        end

        actual = - mean(gradient(log10(abs(difference)), log10(tN)));
        hold on;
        plot(log10(tN), log10(abs(difference)));

        verifyEqual(testCase, actual, expected{n}, 'AbsTol', 0.1);
    end

end

function testOdefunSolveJacobian(testCase)
    explicitOdefun = @(t,y) y;
    t = linspace(0,10,100)';
    y0 = 0.1;
    options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
        optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true)));

    solverList = {@bdf1si, @bdf2si, @bdf3si};
    expectedAccuracy = {1e-8, 1e-8, 1e-8};

    A = y0 / (1 - y0);
    expected = A*exp(t')./(1 + A*exp(t'));
    plot(expected); hold on;
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        solution = solverList{n}(explicitOdefun, @implicitOdefun, t, y0, options);
        actual = solution;
        plot(actual);

        verifyEqual(testCase, actual, expected, 'AbsTol', expectedAccuracy{n});
    end

    function [f, J] = implicitOdefun(t, y)
        f = - y.^2;
        J = - 2 * y;
    end
end

function testConvergenceRatesJacobian(testCase)
    explicitOdefun = @(t,y) y;
    tN = round(logspace(2,2.5,6)); %[100, 110, 120];
    tL = 1;
    y0 = 0.5;
    options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
        optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true)));

    solverList = {@bdf1si, @bdf2si, @bdf3si};
    expected = {2, 3, 4};

    difference = zeros(length(tN),1);
    trueValue = exp(tL)./(1 + exp(tL));
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));

        for m = 1:length(tN)
            t = linspace(0,tL,tN(m))';
            solution = solverList{n}(explicitOdefun, @implicitOdefun, t, y0, options);
            difference(m) = solution(end) - trueValue;
        end

        actual = - mean(gradient(log10(abs(difference)), log10(tN)));
        hold on;
        plot(log10(tN), log10(abs(difference)));

        verifyEqual(testCase, actual, expected{n}, 'AbsTol', 0.1);
    end

    function [f, J] = implicitOdefun(t, y)
        f = - y.^2;
        J = - 2 * y;
    end
end
