function [ R ] = quatToRot( quat )
% quatToR takes unit quaternion to rotation matrix
%   [R] = quatToR(quat)
%   quat= [q0,q1,q2,q3] assumed to follow q0 = cos(angle / 2) while
%   [q1,q2,q3] = axis*sin(angle / 2) for angle/axis expression of R

e0 = quat(1);
e1 = quat(2);
e2 = quat(3);
e3 = quat(4);

R = [1-2*(e2^2+e3^2)   2*(e1*e2-e0*e3)     2*(e1*e3 + e0*e2) ;
     2*(e1*e2+e0*e3)   1-2*(e1^2+e3^2)     2*(e2*e3 - e0*e1) ;
     2*(e1*e3 -e0*e2)  2*(e2*e3+e0*e1)     1-2*(e1^2+e2^2)];
 
end

