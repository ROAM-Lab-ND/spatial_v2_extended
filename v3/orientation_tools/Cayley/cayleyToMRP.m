function m = cayleyToMRP(c,n)
if nargin == 1
    n = 1;
end

if n ==1
    m = c/(1+sqrt(1+c'*c));
else
    m = rotToMRP( calyeyToRot(c,n) );
end

end

