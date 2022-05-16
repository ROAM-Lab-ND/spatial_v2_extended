function q = cayleyToQuat(c,n)

if nargin == 1
    n = 1;
end

if n == 1
    q = 1/sqrt(1+c'*c)*[1 ;c];
else
    q = rotToQuat( cayleyToRot(c,n) );
end

end

