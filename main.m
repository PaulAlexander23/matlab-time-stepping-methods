%MAIN Script to dummy run the various time stepping methods

y0 = 1;
eta = 1;

t = linspace(0,1,10);

solver = 'ab3';
problem = @(t,y) ode(y,eta);

tic;
y = feval(solver,problem,t,y0);
tc = toc;

if isstruct(y)
    t = y.x;
    y = y.y;
end

yexact = y0 * exp(eta * t);
ye = y - yexact;
err =  norm(ye(end,:),'inf');

fprintf('Problem: %s, Solver: %s, Error: %g at t = %gs, Time: %gs\n',func2str(problem),solver,err,t(end),tc);

plot(t,y,t,yexact)

function f = ode(y,eta)
    f = eta * y;
end