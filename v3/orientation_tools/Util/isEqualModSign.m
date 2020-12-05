function flag = isEqualModSign(a,b,tol)

if nargin == 2
    tol = sqrt(eps)/2;
end

flag = isEqual(a,b,tol) | isEqual(a,-b,tol);
end

