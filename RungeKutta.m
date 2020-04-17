classdef RungeKutta < timeStepper
    methods
        function y = run(obj, odefun, t, y0)
            y = rk4(odefun, t, y0);
        end
    end
end