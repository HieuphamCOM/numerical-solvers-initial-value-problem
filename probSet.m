classdef probSet
   properties ( Access = public )
       t0;
       tf;
       dt;
       dtTest;
       t;
       y0;
       yprime;
       yexact;
       

   end
   
   methods
       function prob = probSet(type)
           switch type
               case 1
                   prob.t0 = 0;
                   prob.tf = 5;
                   prob.dt = 0.01;
                   prob.t = prob.t0:prob.dt:prob.tf;
                   prob.y0 = 0;
                   prob.yprime = @(y,t) y - 0.5*exp(t/2)*sin(5*t) + 5*exp(t/2)*cos(5*t);
                   prob.yexact = @(t) exp(t/2).*sin(5*t);
                   prob.dtTest = [0.5 0.1 0.05 0.01 0.005 0.001]';
           end
       end
   end
end