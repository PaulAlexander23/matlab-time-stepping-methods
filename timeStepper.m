classdef timeStepper
    properties
        options
    end
    methods
        function obj = timeStepper()
            obj.options = odeset();
        end
        
        function obj = setAbsTol(obj, absTol)
            obj.options.AbsTol = absTol;
        end
        
        function obj = setRelTol(obj, relTol)
            obj.options.RelTol = relTol;
        end
        
        function obj = setJacobian(obj, jacobian)
            obj.options.Jacobian = jacobian;
        end
    end
    methods (Abstract)
        run(odefun, t, y0)
    end
end