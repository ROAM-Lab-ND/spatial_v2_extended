function q = cayleyToQuat(c)
q = 1/sqrt(1+c'*c)*[1 ;c];
end

