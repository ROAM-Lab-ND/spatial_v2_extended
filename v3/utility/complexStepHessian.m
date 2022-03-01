function [ ddtau_dq, ddtau_dqd,steps] = complexStepHessian(f, x, step)
    if nargin == 2
        step = eps;
    end
    n = length(f(x));
    m = length(x);
    
    H = zeros(n,m,m);
    for ind = 1:m
       e = zeros(m,1);
       e(ind) = 1;
       steps(:,ind) = x+1i*e*step;
       [temp1, temp2 ]= f(x+1i*e*step);
       
        ddtau_dq(:,:,ind) = imag( temp1)/step;  
        ddtau_dqd(:,:,ind)= imag( temp2)/step;  

    end
end