function tests = testImplicitSolvers()
    tests = functiontests(localfunctions);
end

function testExponentialDecayInputDefaults(testCase)
    epsilon = 1;
    odefun = @(t, y) - epsilon * y;
    t = linspace(0,1,10)';
    y0 = 1;
    options = odeset();

    solverList = {@am1, @am2, @am3, @abm1, @abm2, @abm3, @bdf1, @bdf2, @bdf3};
    expectedAccuracy = {2e-2, 4e-4, 3e-5, 2e-2, 4e-4, 7e-5, 2e-2, 3e-2, 3e-2};

    expected = exp(- epsilon * t(end));
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        solution = solverList{n}(odefun, t, y0, options);
        actual = solution.y(end);

        verifyEqual(testCase, actual, expected, 'AbsTol', expectedAccuracy{n});
    end
end

function testExponentialDecayOneOutput(testCase)
    epsilon = 1;
    odefun = @(t, y) - epsilon * y;
    t = linspace(0,1,10)';
    y0 = 1;
    options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
        optimoptions('fsolve', 'Display', 'off')));

    solverList = {@am1, @am2, @am3, @abm1, @abm2, @abm3, @bdf1, @bdf2, @bdf3};
    expectedAccuracy = {2e-2, 4e-4, 3e-5, 2e-2, 4e-4, 7e-5, 2e-2, 3e-2, 3e-2};

    expected = exp(- epsilon * t(end));
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        solution = solverList{n}(odefun, t, y0, options);
        actual = solution.y(end);

        verifyEqual(testCase, actual, expected, 'AbsTol', expectedAccuracy{n});
    end
end

function testExponentialDecayTwoOutputs(testCase)
    epsilon = 1;
    odefun = @(t, y) - epsilon * y;
    t = linspace(0,1,10)';
    y0 = 1;
    options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
        optimoptions('fsolve', 'Display', 'off')));

    solverList = {@am1, @am2, @am3, @abm1, @abm2, @abm3, @bdf1, @bdf2, @bdf3};
    expectedAccuracy = {2e-2, 4e-4, 3e-5, 2e-2, 4e-4, 7e-5, 2e-2, 3e-2, 3e-2};

    expected = exp(- epsilon * t(end));
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        [t, y] = solverList{n}(odefun, t, y0, options);
        actual = y(end);

        verifyEqual(testCase, actual, expected, 'AbsTol', expectedAccuracy{n});
    end
end

function testConvergenceRates(testCase)
    epsilon = 1;
    odefun = @(t, y) - epsilon * y;
    tN = round(logspace(2,2.5,6)); %[100, 110, 120];
    tL = 1;
    y0 = 1;
    options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
        optimoptions('fsolve', 'Display', 'off')));

    solverList = {@am1, @am2, @am3, @abm1, @abm2, @abm3, @bdf1, @bdf2, @bdf3};
    expected = {1, 2, 3, 1, 2, 3, 1, 2, 3};

    difference = zeros(length(tN),1);
    trueValue = exp(- epsilon * tL);
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));

        for m = 1:length(tN)
            t = linspace(0,tL,tN(m))';
            solution = solverList{n}(odefun, t, y0, options);
            difference(m) = solution.y(end) - trueValue;
        end

        actual = - mean(gradient(log10(abs(difference)), log10(tN)));
        % hold on;
        % plot(log10(tN), log10(abs(difference)));

        verifyEqual(testCase, actual, expected{n}, 'AbsTol', 0.1);
    end
end

function testExponentialDecayJacobian(testCase)
    epsilon = 1;
    odefun = @(t, y) - epsilon * y;
    odejac = @(t, y) - epsilon;
    t = linspace(0,1,10)';
    y0 = 1;
    optimmethod = @(fun, x0) fsolve(fun, x0, ...
            optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient',true));
    options = struct( ...
        'optimmethod', optimmethod ...
        ,'Jacobian', odejac ...
        );

    solverList = {@bdf1};
    expectedAccuracy = {2e-2};

    expected = exp(- epsilon * t(end));
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        solution = solverList{n}(odefun, t, y0, options);
        actual = solution.y(end);

        verifyEqual(testCase, actual, expected, 'AbsTol', expectedAccuracy{n});
    end
end
