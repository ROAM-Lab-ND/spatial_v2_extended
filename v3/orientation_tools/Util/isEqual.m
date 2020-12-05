function flag = isEqual(a,b,tol)

if nargin == 2
    tol = sqrt(eps)/2;
end
d = a-b;
flag = max(abs(d(:)))<=tol;

end

