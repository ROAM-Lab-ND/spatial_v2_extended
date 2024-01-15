function mat = icrf(f)

mat = [-skew2(f(1:3)) -skew2(f(4:6));-skew2(f(4:6)) zeros(3)];
 
