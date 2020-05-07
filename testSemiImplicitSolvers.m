function tests = testSemiImplicitSolvers()
    tests = functiontests(localfunctions);
end

function testOutputStruct(testCase)
    explicitOdefun = @(t,y) 1;
    implicitOdefun = @(t,y) -1;
    odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
    t = linspace(0,1,10)';
    y0 = 1;
    myoptimoptions = optimoptions('fsolve', 'Display', 'off');
    options = struct(...
        'optimmethod', @fsolve, ...
        'optimoptions', myoptimoptions);

    solverList = {@ab1be, @ab2be, @ab3cn, @bdf1si, @bdf2si, @bdf3si, @bdf4si};

    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        solution = solverList{n}(odefun, t, y0, options);

        expected = t;
        actual = solution.x;
        verifyEqual(testCase, actual, expected);

        actual = solution.y';
        verifySize(testCase, actual, [length(y0),length(t)]);
    end
end

function testOutputMultiple(testCase)
    explicitOdefun = @(t,y) 1;
    implicitOdefun = @(t,y) -1;
    odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
    t = linspace(0,1,10)';
    y0 = 1;
    myoptimoptions = optimoptions('fsolve', 'Display', 'off');
    options = struct(...
        'optimmethod', @fsolve, ...
        'optimoptions', myoptimoptions);

    solverList = {@ab1be, @ab2be, @ab3cn, @bdf1si, @bdf2si, @bdf3si, @bdf4si};

    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        [actualT, actualY] = solverList{n}(odefun, t, y0, options);

        expected = t;
        verifyEqual(testCase, actualT, expected);

        verifySize(testCase, actualY', [length(y0),length(t)]);
    end
end

function testOdefunSolveError(testCase)
    explicitOdefun = @(t,y) y;
    implicitOdefun = @(t,y) -y^2;
    odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
    t = linspace(0,7,100)';
    y0 = 0.01;
    myoptimoptions = optimoptions('fsolve', 'Display', 'off');
    options = struct(...
        'optimmethod', @fsolve, ...
        'optimoptions', myoptimoptions);

    solverList = {@ab1be, @ab2be, @ab3cn, @bdf1si, @bdf2si, @bdf3si, @bdf4si};
    expectedAccuracy = {6e-2, 3e-2, 4e-2, 6e-2, 4e-3, 9e-4, 1e-8};

    A = y0 / (1 - y0);
    expected = A*exp(t')./(1 + A*exp(t'));
    % plot(expected); hold on;
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        solution = solverList{n}(odefun, t, y0, options);
        actual = solution.y';
        % plot(actual);

        verifyEqual(testCase, actual, expected, 'AbsTol', expectedAccuracy{n});
    end
end

function testOdefunSolveConvergenceRates(testCase)
    explicitOdefun = @(t,y) y;
    implicitOdefun = @(t,y) -y^2;
    odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
    tN = round(logspace(2,2.7,6)); %[100, 110, 120];
    tL = 7;
    y0 = 0.01;
    myoptimoptions = optimoptions('fsolve', 'Display', 'off');
    options = struct(...
        'optimmethod', @fsolve, ...
        'optimoptions', myoptimoptions);

    solverList = {@ab1be, @ab2be, @ab3cn, @bdf1si, @bdf2si, @bdf3si, @bdf4si};
    expected = {1, 1, 2, 1, 2, 3, 4};

    difference = zeros(length(tN),1);
    A = y0 / (1 - y0);
    trueValue = A*exp(tL)./(1 + A*exp(tL));
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

        verifyEqual(testCase, actual, expected{n}, 'AbsTol', 0.2);
    end

end

function testOdefunSolveJacobian(testCase)
    explicitOdefun = @(t,y) y;
    implicitOdefun = @(t,y) -y.^2;
    odejac = @(t,y) -2*y;
    odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
    t = linspace(0,7,100)';
    y0 = 0.01;
    optimmethod = @fsolve;
    myoptimoptions = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true);
    options = struct(...
        'optimmethod', optimmethod, ...
        'optimoptions', myoptimoptions, ...
        'Jacobian', odejac);

    solverList = {@bdf1si, @bdf2si, @bdf3si};
    expectedAccuracy = {6e-2, 4e-3, 9e-4};

    A = y0 / (1 - y0);
    expected = A*exp(t')./(1 + A*exp(t'));
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        solution = solverList{n}(odefun, t, y0, options);

        actual = solution.y';
        verifyEqual(testCase, actual, expected, 'AbsTol', expectedAccuracy{n});
    end
end

function testDampedOscillatorError(testCase)
    b = -0.1;
    explicitOdefun = @(t,y) [0; -2 * b * y(2)];
    implicitOdefun = @(t,y) [y(2); -y(1)];
    odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
    t = linspace(0,(2.25)*pi,500)';
    y0 = [1;-b];
    myoptimoptions = optimoptions('fsolve', 'Display', 'off');
    options = struct(...
        'optimmethod', @fsolve, ...
        'optimoptions', myoptimoptions);

    solverList = {@ab1be, @ab2be, @ab3cn, @bdf1si, @bdf2si, @bdf3si};
    expectedAccuracy = {9e-2, 9e-2, 2e-1, 9e-2, 2e-3, 3e-4};

    B = sqrt(1-b^2);
    expected = [exp(-b * t') .* cos(B*t'); -b * exp(-b * t') .* cos(B*t') - B*exp(-b * t') .* sin(B*t')];
    % figure(1);
    % plot(t, expected(1,:));
    % figure(2);
    % plot(t, expected(2,:));
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));
        solution = solverList{n}(odefun, t, y0, options);
        actual = solution.y';

        % figure(1);
        % hold('on'); plot(t, actual(1,:)); hold('off');
        % figure(2);
        % hold('on'); plot(t, actual(2,:)); hold('off');

        verifyEqual(testCase, actual, expected, 'AbsTol', expectedAccuracy{n});
    end
end

function testDampedOscillatorConvergence(testCase)
    b = -0.1;
    explicitOdefun = @(t,y) [0; -2 * b * y(2)];
    implicitOdefun = @(t,y) [y(2); -y(1)];
    odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
    tN = round(logspace(2.5,3.5,4)); %[100, 110, 120];
    tL = (2.25)*pi;
    y0 = [1;-b];
    myoptimoptions = optimoptions('fsolve', 'Display', 'off');
    options = struct(...
        'optimmethod', @fsolve, ...
        'optimoptions', myoptimoptions);

    solverList = {@ab1be, @ab2be, @ab3cn, @bdf1si, @bdf2si, @bdf3si};
    expected = {1, 1, 2, 1, 2, 3};

    difference = zeros(length(y0), length(tN));
    B = sqrt(1-b^2);
    trueValue = [exp(-b * tL) .* cos(B*tL); -b * exp(-b * tL) .* cos(B*tL) - B*exp(-b * tL) .* sin(B*tL)];
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));

        for m = 1:length(tN)
            t = linspace(0,tL,tN(m))';
            solution = solverList{n}(odefun, t, y0, options);
            difference(:,m) = solution.y(end,:)' - trueValue;
        end

        actual = - mean(gradient(log10(abs(difference)), log10(tN),1),2);
        % figure(1);
        % hold on;
        % plot(log10(tN), log10(abs(difference(1,:))));
        % hold off;
        % figure(2);
        % hold on;
        % plot(log10(tN), log10(abs(difference(2,:))));
        % hold off;

        verifyEqual(testCase, actual, [expected{n}; expected{n}], 'AbsTol', 0.1);
    end

end
