function L = quatL(quat)
% quatL takes a quaternion q1 and returns a 4x4 matrix such that
%   quatL(q1)*q2 = q1*q2 (where * means quaternion multiplication)


sca = quat(1);
vec = quat(2:4);
L = [sca -vec.' ; vec sca*eye(3)+skew(vec)];

end

