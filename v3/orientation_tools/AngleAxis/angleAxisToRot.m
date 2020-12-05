function R = angleAxisToRot(angle_axis)
%R = expm(skew(angle_axis));

sinc2 = @(z) sinc(z/pi);
e = angle_axis;
phi = sqrt(e'*e);
I = eye(3);

% Rodrigues formula, without having to sorry about divide by zero.
R = cos(phi)*I + sinc2(phi/2)^2*e*e'/2 + sinc2(phi)*skew(e);

end

