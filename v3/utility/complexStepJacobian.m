function J = complexStepJacobian(f, x, step)
    if nargin == 2
        step = eps;
    end
    n = length(f(x));
    m = length(x);
    
    J = zeros(n,m);
    for ind = 1:m
       e = zeros(m,1);
       e(ind) = 1;
       J(:,ind) = imag( f(x+i*e*step))/step;  
    end
end