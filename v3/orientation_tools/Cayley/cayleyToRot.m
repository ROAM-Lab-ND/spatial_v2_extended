function R = cayleyToRot(c , n)

% x = c(1);
% y = c(2);
% z = c(3);
% % See https://en.wikipedia.org/wiki/Rotation_matrix#Skew_parameters_via_Cayley's_formula
% R = 1/(1+x^2+y^2+z^2)*[1+x^2-y^2-z^2 2*x*y-2*z 2*y+2*x*z; ... 
% 2*x*y+2*z 1-x^2+y^2-z^2 2*y*z-2*x ; ... 
% 2*x*z-2*y 2*x+2*y*z 1-x^2-y^2+z^2];

if nargin == 1
    n = 1;
end

if n == 1
    cc = c'*c;
    R = ( eye(3)*(1-cc)+2*c*c' + 2*skew(c) ) / (1+cc);
else
    R = (eye(3)-skew(-c))^n/(eye(3) + skew(-c))^n;
end

end

