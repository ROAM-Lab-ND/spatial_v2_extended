function L = quatL(quat)

sca = quat(1);
vec = quat(2:4);
L = [sca -vec' ; vec sca*eye(3)+skew(vec)];

end

