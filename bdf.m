function y = bdf(odefun,t,y0)
    n = length(t);
    y = zeros([size(y0),n]);
    
    y(:,:,1) = y0;
    h = y0;
    hm1 = y0;
    hm2 = y0;
    
    for i = 2:n
        dt = t(i)-t(i-1);
        
        [h, hm1,  hm2] = step_newton(odefun,t(i),dt,h,hm1,hm2);
        
        y(:,:,i) = h;
    end
    
    y = squeeze(y);
    
    function [h, hm1,  hm2] = step_newton(odefun,t,dt,h,hm1,hm2)
        
        temp = h;
        
        error = 1;
        m = 1;
        
        while error > 1e-7 && m <= 20
            res =  h - bdf2(odefun,t,dt,h,hm1,hm2);
            error = norm(res);
            jac = jacobian(odefun,t,dt,h);
            h = h-jac\res;
            m = m + 1;
        end
        
        hm2 = hm1;
        hm1 = temp;
    end
    
    function out = bdf1(odefun,t,dt,h,hm1) %DT MUST BE FIXED
        out = hm1 + dt * odefun(t, h);
    end
    
    function out = bdf2(odefun,t,dt,h,hm1,hm2) %DT MUST BE FIXED
        out = 4/3 * hm1 - 1/3 * hm2 + 2/3 * dt * odefun(t, h);
    end
    
    function out = bdf3(odefun,t,dt,h,hm1,hm2,hm3) %DT MUST BE FIXED
        out = 18/11 * hm1 - 9/11 * hm2 + 2/11 * hm3 + 6/11 * dt * odefun(t, h);
    end
    
    function out = jacobian(odefun,t,dt,h)
        %J_control = [I/(2*delta_t/3)+Dx*q_h];
    end
end