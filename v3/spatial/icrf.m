function mat = icrf(f)

mat = [-skew(f(1:3)) -skew(f(4:6));-skew(f(4:6)) zeros(3)];
 
