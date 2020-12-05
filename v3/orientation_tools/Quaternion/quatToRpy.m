function [ rpy ] = quatToRpy( quat )
% quatToRpy takes unit quaternion to earth-fixed sequenial roll-pitch-yaw
%           euler angles
%   [rpy] = quatToRpy(quat)
%   quat= [q0,q1,q2,q3] assumed to follow q0 = cos(angle / 2) while
%   [q1,q2,q3] = axis*sin(angle / 2) for angle/axis expression of R
%   Note: earth-fixed roll-pitch-yaw is the same as body-fixed
%   yaw-pitch-roll sequence 

q0 = quat(1);
q1 = quat(2);
q2 = quat(3);
q3 = quat(4);

rpy = [atan2(2*q2*q3+2*q0*q1,q3^2-q2^2-q1^2+q0^2);
       -asin(2*q1*q3-2*q0*q2);
       atan2(2*q1*q2+2*q0*q3,q1^2+q0^2-q3^2-q2^2)];
end

