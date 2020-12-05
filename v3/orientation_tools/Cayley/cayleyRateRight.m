function cDot = cayleyRateRight(c,w)
I = eye(3);    
cDot = 1/2*( I + skew(c) + c*c')*w;
end

