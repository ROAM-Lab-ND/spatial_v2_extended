function J = finiteDiffJacobian(f, x, step)
    if nargin == 2
        step = sqrt(eps);
    end
    f0 = f(x);
    n = length(f(x));
    m = length(x);
    
    J = zeros(n,m);
    for ind = 1:m
       e = zeros(m,1);
       e(ind) = 1;
       J(:,ind) = real(f(x+e*step) - f0)/step;  
    end
end