function [J, evals, steps] = complexStepJacobian(f, x, step)
    if nargin == 2
        step = eps;
    end
    n = length(f(x));
    m = length(x);
    
    J = zeros(n,m);
    for ind = 1:m
       e = zeros(m,1);
       e(ind) = 1;
       steps(:,ind) = x+1i*e*step;
       evals(:,ind) = f(x+1i*e*step);
       
       J(:,ind) = imag( f(x+1i*e*step))/step;  
    end
end