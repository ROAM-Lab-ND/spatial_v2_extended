function [J, evals, steps] = complexStepJacobian(f, x, step)
    if nargin == 2
        step = eps;
    end
    f0 = f(x);
    f0_vec = f0(:);
    
    n = length(f0_vec);
    m = length(x);
    
    J = zeros(n,m);
    for ind = 1:m
       e = zeros(m,1);
       e(ind) = 1;
       steps(:,ind) = x+1i*e*step;
       evals(:,ind) = reshape( f(x+1i*e*step), [], 1);
       
       J(:,ind) = imag( evals(:,ind) )/step;  
    end
    
    evals = reshape(evals, [size(f0) length(x)]);
    J     = reshape(J, [size(f0) length(x)]);
	if size(J,2) == 1
        J = squeeze(J);
	end
end