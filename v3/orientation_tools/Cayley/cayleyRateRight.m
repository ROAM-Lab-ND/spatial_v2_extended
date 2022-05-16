function cDot = cayleyRateRight(c,w,n)
    if nargin == 2
        n = 1;
    end
    I = eye(3);    

    if n ==1 
        cDot = 1/2*( I + skew(c) + c*c')*w;
    else
        cDot = getRateRight(c, w,n);
    end

end

%%%%%%%%%%%%%%%%
function [dtau, A] = getRateRight(tau,w, n)
    for i = 1:3
        ei = zeros(3,1);
        ei(i) = 1;
        A(:,i) = getOmegaRight(tau,ei,n);
    end
    dtau = A\w;
end

function w = getOmegaRight(tau,dtau,n)
    [w, R] = getOmegaLeft(tau, dtau, n );
    w = R'*w;
end

function [w, R] = getOmegaLeft(tau, dtau, n )
    T  = -skew(tau);
    dT = -skew(dtau);
    [sk_w,  R] = derivHelper(T,dT,n);
    w = skew(sk_w);
end

function [sk_w,R] = derivHelper(T,dT,n)
   A = ( eye(3) - T )^n;
   B = ( eye(3) + T )^n;
   R = A/B;
   
   D1 = differ(eye(3)-T,n, -dT);
   D2 = R*differ(eye(3)+T,n,dT);
   sk_w = (D1 - D2)/A;
end

function D = differ(A,n, dA)
    D = 0*A;
    for j = 0:(n-1)
        D = D + A^j * dA * A^(n-j-1);
    end
end
