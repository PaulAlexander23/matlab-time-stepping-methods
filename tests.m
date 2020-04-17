function tests = tests()
    tests = functiontests(localfunctions);
end

function testAdamsBashford2BackwardEuler(testCase)
    t = linspace(0, 10, 100);
    y0 = 1;
    odefunLinear = @(t, y) 2 * y;
    odefunNonlinear = @(t, y) - y.^2 + exp(- 4 * t);
    options = optimoptions('fsolve','Display', 'off');
    optimmethod = @(fun, x0) fsolve(fun, x0, options);
    
    actual = ab2be(odefunLinear, odefunNonlinear, t, y0, optimmethod);
    expected = 2 - exp(-2 * t);
    
    verifyEqual(testCase, actual, expected, 'RelTol', 10*(t(2)-t(1))^2)
end

function testAdamsBashford3CrankNicholson(testCase)
    t = linspace(0, 10, 100);
    y0 = 1;
    odefunLinear = @(t) 2;
    odefunNonlinear = @(t, y) - y.^2 + exp(- 4 * t);
    options = struct("Linear", odefunLinear);
    
    actual = ab3cn(odefunNonlinear, t, y0, options);
    expected = 2 - exp(-2 * t);
    
    verifyEqual(testCase, actual, expected, 'RelTol', 6e-2)
end

function testAdamsBashford3CrankNicholsonOrder(testCase)
    tN = 2.^(3:10);
    N = length(tN);
    maxError = zeros(1,N);
    
    for k = 1:N
        t = linspace(0, 10, tN(k));
        y0 = 1;
        odefunLinear = @(t) 2;
        odefunNonlinear = @(t, y) - y.^2 + exp(- 4 * t);
        options = struct("Linear", odefunLinear);

        actual = ab3cn(odefunNonlinear, t, y0, options);
        expected = 2 - exp(-2 * t);
        error = actual - expected;
        maxError(k) = max(abs(error));
    end
    
    expectedConvergence = 2;
    
    actual = log2(maxError);
    expected = -expectedConvergence * (log2(tN) - log2(tN(end))) + log2(maxError(end));
    
    cutoff = 6;
    verifyEqual(testCase, actual(cutoff:end), expected(cutoff:end), 'RelTol', 1e-2) 
end
