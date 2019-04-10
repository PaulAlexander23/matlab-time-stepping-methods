%MAIN Script to dummy run the various time stepping methods

y0 = 1;
eta = 1;

t = linspace(0,1,10);

solver = 'bdf';
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

%% Problem with a jacobian

y0 = [2,3]';

eta = 1;
t = linspace(0,10,101);

% options = optimoptions('fsolve','SpecifyObjectiveGradient',true);
% solver = @(odefun,t,y0) bdf3(odefun,t,y0,@(fun,x0) fsolve(fun,x0,options));
% solver = @(odefun,t,y0) bdf1(odefun,t,y0,@newton);
% options = odeset('Jacobian',@(t,y) jac2(y,eta));
% solver = @(problem,t,y0) ode15s(problem,t,y0,options);
% solver = @ode15s;
options = struct(...'Jacobian',@(t,y) jac2(y,eta), ...
    'ExplicitFcn', @expf, ...
    'Order', 3);
solver = 'bdf';
problem = @(t,y) fun2(y,eta);

tic
y = feval(solver,problem,t,y0,options);
toc

if isstruct(y)
    t = y.x;
    y = y.y;
end

plot(t,y)

%% Problem with an explicit function

y0 = [2,3]';
t = linspace(0,10,101);

optimmethod = @fsolve;
tic
y = bdf1_semi_implicit(@(t,y) fun2(y,1), t, y0, optimmethod, @expf);
toc

plot(t,y)

function [f,J] = fun(y)
    f = [(y(1)-1).^2 + (y(1)-1).^2;(y(1)-1).*(y(2)-1)];
    if nargout > 1
        J = [2*(y(1)-1),2*(y(2)-1);y(2)-1,y(1)-1];
    end
end

function [f, j] = fun2(y,eta)
    f = -[y(1)-1 + eta*(y(1)-1).*(y(2)-1);...
        y(2)-1 - eta*(y(1)-1).*(y(2)-1)];
%     if nargout > 1
%         j = -[1 + eta*(y(2)-1), eta*(y(1)-1);
%             - eta*(y(2)-1), 1 + eta*(y(1)-1)];
%     end
end

function j = jac2(y,eta)
    j = -[1 + eta*(y(2)-1), eta*(y(1)-1);
        - eta*(y(2)-1), 1 + eta*(y(1)-1)];
end
function f = ode(y,eta)
    f = eta * y;
end

function f = expf(t,y)
    f = -y.^2 + (t - floor(t)).*ones(size(y));
end