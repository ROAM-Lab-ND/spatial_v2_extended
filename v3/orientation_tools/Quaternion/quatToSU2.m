function S = quatToSU2(q)
% Converts a quaternion to a special unitary 2x2 matrix

a = q(1);
b = q(2);
c = q(3);
d = q(4);

i = 1i;

S = [a+i*d -b-i*c ; b-i*c a-i*d];

