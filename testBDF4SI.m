function tests = testBDF4SI()
    tests = functiontests(localfunctions);
end

function testDampedOscillatorConvergencePerfectWindupMaxStep(testCase)
    b = -0.1;
    explicitOdefun = @(t,y) [0; -2 * b * y(2)];
    implicitOdefun = @(t,y) [y(2); -y(1)];
    odefun = struct('explicit', explicitOdefun, 'implicit', implicitOdefun);
    tN = 21;
    tL = (2.25)*pi;
    MaxStep = [5e-3, 1e-3] * tL;
    y0 = [1;-b];
    myoptimoptions = optimoptions('fsolve', 'Display', 'off');
    options = odeset();
    options.optimmethod = @fsolve;
    options.optimoptions = myoptimoptions;

    solverList = {@bdf4si};
    expected = {4};

    difference = zeros(length(y0), length(tN));
    B = sqrt(1-b^2);
    analyticalSolution = @(t) [exp(-b * t) .* cos(B*t); -b * exp(-b * t) .* cos(B*t) - B*exp(-b * t) .* sin(B*t)];
    trueValue = [exp(-b * tL) .* cos(B*tL); -b * exp(-b * tL) .* cos(B*tL) - B*exp(-b * tL) .* sin(B*tL)];
    for n = 1:length(solverList)
        fprintf("Solver: %s,\n", func2str(solverList{n}));

        for m = 1:length(MaxStep)
            options.MaxStep = MaxStep(m);
            t = linspace(0,tL,tN)';
            y0 = analyticalSolution(0:options.MaxStep:options.MaxStep*3);
            % y0 = analyticalSolution(0); Imperfect windup
            solution = solverList{n}(odefun, t, y0, options);
            difference(:,m) = solution.y(end,:)' - trueValue;
        end

        actual = mean(gradient(log10(abs(difference)), log10(MaxStep),1),2);

        % figure(1);
        % hold on;
        % plot(log10(MaxStep), log10(abs(difference(1,:))));
        % hold off;
        % figure(2);
        % hold on;
        % plot(log10(MaxStep), log10(abs(difference(2,:))));
        % hold off;

        verifyEqual(testCase, actual, [expected{n}; expected{n}], 'AbsTol', 0.7);
    end
end

