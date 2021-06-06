function R = quatR(quat)
% quatR takes a quaternion q2 and returns a 4x4 matrix such that
%   quatR(q2)*q1 = q1*q2 (where * means quaternion multiplication)

sca = quat(1);
vec = quat(2:4);
R = [sca -vec.' ; vec sca*eye(3)-skew(vec)];

end

