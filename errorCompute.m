classdef errorCompute
   properties (Access = public)
      errorFirstExplicit
      errorFirstImplicit
      errorSecondExplicit
      errorSecondImplicit
      errorRK4
   end
   
   methods
       function err = errorCompute(prob, soln)
          
           err.errorFirstExplicit = zeros(length(prob.dtTest),1);
           err.errorFirstImplicit = zeros(length(prob.dtTest),1);
           err.errorSecondExplicit = zeros(length(prob.dtTest),1);
           err.errorSecondImplicit = zeros(length(prob.dtTest),1);
           err.errorRK4 = zeros(length(prob.dtTest));     
           for i = 1 : length(prob.dtTest)
               prob.t = prob.t0:prob.dtTest(i):prob.tf;
               % Exact
               soln = soln.computeExact(prob); 
               % First Explicit Forward Euler
               soln = soln.computeFirstExplicit(prob);
               % First Implicit Backward Euler
               soln = soln.computeFirstImplicit(prob);
               % Second Explicit 
               soln = soln.computeSecondExplicit(prob);
               % Second Implicit 
               soln = soln.computeSecondImplicit(prob);
               % R-K-4
               soln = soln.computeRK4(prob);
               
               err.errorFirstExplicit(i,1) = abs(soln.solExact(end)-soln.solFristExplicit(end));
               err.errorFirstImplicit(i,1) = abs(soln.solExact(end)-soln.solFristImplicit(end));
               err.errorSecondExplicit(i,1) = abs(soln.solExact(end)-soln.solSecondExplicit(end));
               err.errorSecondImplicit(i,1) = abs(soln.solExact(end)-soln.solSecondImplicit(end));
               err.errorRK4(i,1) = abs(soln.solExact(end)-soln.solRK4(end));
               
           end
           figure(2)
           loglog(prob.dtTest,err.errorFirstExplicit,...
               prob.dtTest,err.errorFirstImplicit,...
               prob.dtTest,err.errorSecondExplicit,...
               prob.dtTest,err.errorSecondImplicit,...
               prob.dtTest,err.errorRK4)
           legend('ForwardEuler','BackwardEuler','MidpointExplicit','MidpointImplicit','Explicit Runge-Kutta 4th order')
       end
   end
end