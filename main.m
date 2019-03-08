
y0 = 1;
eta = 1;

t = linspace(0,1,10)';

f = @(t, y) eta * y;

yexact = y0 * exp(eta * t);

y = abm2(f,t,y0);

err = y - yexact;

fprintf("Error at t = %.1f: %g\n", t(end), norm(err(end,:)));

plot(t,y,t,yexact)