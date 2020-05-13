function tests = testExplicitSolvers()
    tests = functiontests(localfunctions);
end

function testExponentialDecaySingleOutput(testCase)
    y0 = 1;
    epsilon = 1;
    odefun = @(t, y) - epsilon * y;
    t = linspace(0,1,10)';

    solverList = {@ab1, @ab2, @rk4};
    expectedAccuracy = {1e-1, 1e-2, 1e-6};

    expected = exp(- epsilon * t(end));
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        solution = solverList{n}(odefun, t, y0);
        actual = solution.y(end)';

        verifyEqual(testCase, actual, expected, 'AbsTol', expectedAccuracy{n});
    end
end

function testExponentialDecayTwoOutputs(testCase)
    y0 = 1;
    epsilon = 1;
    odefun = @(t, y) - epsilon * y;
    t = linspace(0,1,10)';

    solverList = {@ab1, @ab2, @rk4};
    expectedAccuracy = {1e-1, 1e-2, 1e-6};

    expected = exp(- epsilon * t(end));
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        [t, y] = solverList{n}(odefun, t, y0);
        actual = y(end)';

        verifyEqual(testCase, actual, expected, 'AbsTol', expectedAccuracy{n});
    end
end

function testOutputtingAtSpecifiedTimepointsEqualSteps(testCase)
    y0 = 1;
    epsilon = 1;
    odefun = @(t, y) - epsilon * y;
    t = linspace(0,1,11)';
    options = odeset('MaxStep',0.02);

    solverList = {@ab1, @ab2, @rk4};
    expectedAccuracy = {1e-1, 1e-2, 1e-6};

    expected = exp(- epsilon * t(end));
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        solution = solverList{n}(odefun, t, y0, options);
        actual = solution.y(end)';

        verifyEqual(testCase, actual, expected, 'AbsTol', expectedAccuracy{n});
    end
end

function testOutputtingAtSpecifiedTimepointsUnequalSteps(testCase)
    y0 = 1;
    epsilon = 1;
    odefun = @(t, y) - epsilon * y;
    t = linspace(0,1,11)';
    options = odeset('MaxStep',0.03);

    solverList = {@ab1, @ab2, @rk4};
    expectedAccuracy = {1e-1, 1e-2, 1e-6};
    errors = {false, true, false};

    expected = exp(- epsilon * t(end));
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        try
            solution = solverList{n}(odefun, t, y0, options);
            actual = solution.y(end)';

            verifyEqual(testCase, actual, expected, 'AbsTol', expectedAccuracy{n});
            verifyTrue(testCase, ~errors{n});
        catch testError
            verifyTrue(testCase, errors{n});
        end
    end
end

function testConvergenceRates(testCase)
    y0 = 1;
    epsilon = 1;
    odefun = @(t, y) - epsilon * y;
    tN = round(logspace(2,3,6)); %[100, 110, 120];
    tL = 1;

    solverList = {@ab1, @ab2, @rk4};
    expected = {1, 2, 4};

    difference = zeros(length(tN),1);
    trueValue = exp(- epsilon * tL);
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));


        for m = 1:length(tN)
            t = linspace(0,tL,tN(m))';
            solution = solverList{n}(odefun, t, y0);
            difference(m) = solution.y(end)' - trueValue;
        end

        actual = - mean(gradient(log10(abs(difference)), log10(tN)));
        % hold on;
        % plot(log10(tN), log10(abs(difference)));

        verifyEqual(testCase, actual, expected{n}, 'AbsTol', 0.1);
    end

end
