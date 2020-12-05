function ang_ax_dot = angleAxisRateLeft(ang_ax,w)
I = eye(3);  

sinc2 = @(z) sinc(z/pi);

phi = norm(ang_ax);
ax  = ang_ax / phi;
if phi < sqrt(eps)
    ax = [1 0 0]';
end

% This term is zero at 0, and has zero derivative at 0. 
% So when phi< sqrt(eps), this term is zero to machine precision
% If you were to divide it by phi^2, then you could use ang_ax in place of
% ax below
term = 1-cos(phi/2)/sinc2(phi/2); 

ang_ax_dot = (I - 1/2*skew(ang_ax) + term*skew(ax)^2)*w;

end

