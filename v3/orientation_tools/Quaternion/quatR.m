function R = quatR(quat)

sca = quat(1);
vec = quat(2:4);
R = [sca -vec' ; vec sca*eye(3)-skew(vec)];

end

