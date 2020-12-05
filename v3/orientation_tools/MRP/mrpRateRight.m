function mDot = mrpRateRight(m,w)
I = eye(3);    
s = m'*m;
mDot = 1/4*[(1-s)*I+2*skew(m)+2*m*m']*w;
end

