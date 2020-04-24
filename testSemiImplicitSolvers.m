function tests = testSemiImplicitSolvers()
    tests = functiontests(localfunctions);
end

function testOdefunSolveOneOutput(testCase)
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
        actual = solution.y';
        plot(actual);

        verifyEqual(testCase, actual, expected, 'AbsTol', expectedAccuracy{n});
    end
end

function testOdefunSolveTwoOutputs(testCase)
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
        [t, y] = solverList{n}(explicitOdefun, implicitOdefun, t, y0, options);
        actual = y';
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
            difference(m) = solution.y(end) - trueValue;
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
        actual = solution.y';
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
            difference(m) = solution.y(end) - trueValue;
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

function testOdefunSolveOneOutputOscillator(testCase)
    b = -0.1;
    explicitOdefun = @(t,y) [0; -2 * b * y(2)];
    implicitOdefun = @(t,y) [y(2); -y(1)];
    t = linspace(0,(7.25)*pi,500)';
    y0 = [1;-b];
    options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
        optimoptions('fsolve', 'Display', 'off')));

    solverList = {@ab1be, @ab2be, @ab3cn, @bdf1si, @bdf2si, @bdf3si};
    expectedAccuracy = {1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8};

    B = sqrt(1-b^2);
    expected = [exp(-b * t') .* cos(B*t'); -b * exp(-b * t') .* cos(B*t') - B*exp(-b * t') .* sin(B*t')];
    figure(1);
    plot(t, expected(1,:));
    figure(2);
    plot(t, expected(2,:));
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        solution = solverList{n}(explicitOdefun, implicitOdefun, t, y0, options);
        actual = solution.y';

        figure(1);
        hold('on'); plot(t, actual(1,:)); hold('off');
        figure(2);
        hold('on'); plot(t, actual(2,:)); hold('off');

        verifyEqual(testCase, actual, expected, 'AbsTol', expectedAccuracy{n});
    end
end

function testConvergenceRatesOscillator(testCase)
    b = -0.1;
    explicitOdefun = @(t,y) [0; -2 * b * y(2)];
    implicitOdefun = @(t,y) [y(2); -y(1)];
    tN = round(logspace(2.5,3.5,6)); %[100, 110, 120];
    tL = (7.25)*pi;
    y0 = [1;-b];
    options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
        optimoptions('fsolve', 'Display', 'off')));

    solverList = {@ab1be, @ab2be, @ab3cn, @bdf1si, @bdf2si, @bdf3si};
    expected = {[2;2], [2;2], [3;3], [2;2], [3;3], [4;4]};

    difference = zeros(length(y0), length(tN));
    B = sqrt(1-b^2);
    trueValue = [exp(-b * tL) .* cos(B*tL); -b * exp(-b * tL) .* cos(B*tL) - B*exp(-b * tL) .* sin(B*tL)];
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));

        for m = 1:length(tN)
            t = linspace(0,tL,tN(m))';
            solution = solverList{n}(explicitOdefun, implicitOdefun, t, y0, options);
            difference(:,m) = solution.y(end,:)' - trueValue;
        end

        actual = - mean(gradient(log10(abs(difference)), log10(tN),1),2);
        figure(1);
        hold on;
        plot(log10(tN), log10(abs(difference(1,:))));
        hold off;
        figure(2);
        hold on;
        plot(log10(tN), log10(abs(difference(2,:))));
        hold off;

        verifyEqual(testCase, actual, expected{n}, 'AbsTol', 0.1);
    end

end
