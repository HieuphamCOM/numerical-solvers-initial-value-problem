classdef numericalMethods
    
    properties (Access = public)
       solExact; % Exact solution
       solFristExplicit; % First Explicit Forward Euler
       solFristImplicit; % First Implicit Backward Euler
       solSecondExplicit; % Second Exlicit 
       solSecondImplicit; % Second Implicit
       solRK4; % R-K 4
       
    end % end of properties
    
    methods
        function soln = numericalMethods()
            
        end
        function soln = computeExact(soln, prob)
            soln.solExact = zeros(length(prob.t),1);
            for i = 1 : length(prob.t)
                soln.solExact(i) = prob.yexact(prob.t(i));                
            end

        end % end of computeExact
        
        function soln = computeFirstExplicit(soln, prob)
            soln.solFristExplicit = zeros(length(prob.t),1);
            soln.solFristExplicit(1) = prob.y0;
            for i = 2 : length(prob.t) 
               soln.solFristExplicit(i,1) = soln.solFristExplicit(i-1)...
                   + prob.dt*prob.yprime(soln.solFristExplicit(i-1),prob.t(i-1)); 
            end        
        end % end of computeFirstExplicit

        function soln = computeFirstImplicit(soln, prob)
            soln.solFristImplicit = zeros(length(prob.t),1);
            soln.solFristImplicit(1) = prob.y0;
            syms y;
            for i = 2 : length(prob.t) 
               soln.solFristImplicit(i,1) =...
                   solve(y - prob.dt * prob.yprime(y, prob.t(i)) == soln.solFristImplicit(i-1,1), y );
            end 
        end
        
        function soln = computeSecondExplicit(soln, prob)
            soln.solSecondExplicit = zeros(length(prob.t),1);
            soln.solSecondExplicit(1) = prob.y0;
            for i =2 : length(prob.t)
                mid_tempt = soln.solSecondExplicit(i-1)...
                    + (prob.dt/2)* prob.yprime(soln.solSecondExplicit(i-1),prob.t(i-1));
                soln.solSecondExplicit(i) = soln.solSecondExplicit(i-1)...
                    + prob.dt * prob.yprime(mid_tempt, prob.t(i-1) + prob.dt/2);
            end
        end
        
        function soln = computeSecondImplicit(soln, prob)
            soln.solSecondImplicit = zeros(length(prob.t),1);
            soln.solSecondImplicit(1) = prob.y0;
            for i = 2 : length(prob.t)
               y_P_tempt = soln.solSecondImplicit(i-1)...
                   + prob.dt * prob.yprime(soln.solSecondImplicit(i-1),prob.t(i-1)); 
               soln.solSecondImplicit(i) = soln.solSecondImplicit(i-1)...
                   + prob.dt/2 ...
               * ( prob.yprime(soln.solSecondImplicit(i-1),prob.t(i-1)) + prob.yprime(y_P_tempt,prob.t(i)) );
            end
        end
        
        function soln = computeRK4(soln, prob)
            soln.solRK4 = zeros(length(prob.t),1);
            soln.solRK4(1) = prob.y0;
            for i = 2: length(prob.t)
                Dy1 = prob.dt*prob.yprime(soln.solRK4(i-1),prob.t(i-1));
                Dy2 = prob.dt*prob.yprime(soln.solRK4(i-1)+Dy1/2,prob.t(i-1)+prob.dt/2);
                Dy3 = prob.dt*prob.yprime(soln.solRK4(i-1)+Dy2/2,prob.t(i-1)+prob.dt/2);
                Dy4 = prob.dt*prob.yprime(soln.solRK4(i-1)+Dy3  ,prob.t(i-1)+prob.dt );
                
                soln.solRK4(i) = soln.solRK4(i-1) + (1/6) * (Dy1 + 2 * Dy2 + 2 * Dy3 + Dy4);
            end
        end
        
        function soln = solveNumerical(soln,prob)
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
           
           % plotting
           figure(1)
           plot(prob.t, soln.solExact,...
               prob.t, soln.solFristExplicit,...
               prob.t, soln.solFristImplicit,...
               prob.t, soln.solSecondExplicit,...
               prob.t, soln.solSecondImplicit,...
               prob.t, soln.solRK4);
           legend('Exact','First Explicit Euler',...
               'First Implicit Euler','Second Explicit',...
               'Second Implicit','RK4')

        end

    end % end of methods
    
end % end of classdef