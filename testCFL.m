function tests = tests()
    tests = functiontests(localfunctions);
end

function testCourantNumberLinearConvection(testCase)
    dx = 0.01;
    dt = 0.0001;

    x = (dx:dx:1)';
    t = 0:dt:0.1;
    diff = @(y, x) (circshift(y, 1) - circshift(y, -1))/(2*x(1));
    alpha = 0.9;
    odefun = @(t, y) - alpha .* diff(y, x);
    y0 = cos(2*pi*x);

    y = ab1(odefun, t, y0);
    
    hold on;
    f = odefun(t(end-1), y(:, end-1));
    u = f ./ - diff(y(:, end-1), x);
    plot(x, f);
    plot(x, u);
    plot(x, y(:, end-1:end));

    expected = alpha;
    actual = max(abs(u));
    %actual = courant(max(abs(u)), dt, dx)

    verifyEqual(testCase, actual, expected, 'AbsTol', eps);
end

function testCourantNumberLinearConvectionWithDiffusion(testCase)
    dx = 0.01;
    dt = 0.01;

    x = (dx:dx:1)';
    t = 0:dt:0.1;
    diff = @(y, x) (circshift(y, 1) - circshift(y, -1))/(2*x(1));
    alpha = 1;
    odefun = @(t, y) - alpha .* diff(y, x) + diff(diff(y, x), x);
    y0 = cos(2*pi*x);

    [~, y] = ode45(odefun, t, y0);
    y = y';
    
    hold on;
    f = odefun(t(end-1), y(:, end-1));
    u = f ./ - diff(y(:, end-1), x);
    plot(x, f);
    plot(x, u);
    plot(x, y(:, end-1:end));

    expected = alpha;
    actual = max(abs(u));
    %actual = courant(max(abs(u)), dt, dx);

    verifyEqual(testCase, actual, expected);
end

function testCourantNumber(testCase)
    dx = 0.0001;
    dt = 0.1;

    x = (dx:dx:1)';
    t = [0,dt];
    diff = @(y, x) (circshift(y, 1) - circshift(y, -1))/(2*x(1));
    odefun = @(t, y) - y .* diff(y, x);
    y0 = 1.1*cos(2*pi*x);

    y = ab1(odefun, t, y0);

    size(y)
    hold on;
    f = odefun(t(end-1), y(:, end-1));
    u = f ./ -diff(y(:, end-1), x);
    plot(x, f);
    plot(x, u);
    plot(x, y(:,end));

    expected = 0.1;
    actual = max(abs(u));
    %actual = courant(max(abs(u)), dt, dx);

    verifyEqual(testCase, actual, expected);
end
