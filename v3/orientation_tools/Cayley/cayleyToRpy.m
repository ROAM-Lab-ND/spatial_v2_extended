function rpy = cayleyToRpy(c,n)

if nargin == 1
    n  = 1;
end

rpy = rotToRpy( cayleyToRot(c,n) );

end

