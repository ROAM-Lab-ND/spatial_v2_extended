function [ cayley ] = rotToCayley( R )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%     angAx = rotToAngleAxis(R);
%     ang = norm(angAx);
%     if ang > 0
%         ax  = angAx / ang;
%     else
%         ax = [1 0 0]';
%     end
%     cayley = [tan(ang/2) * ax];

I = eye(3);
cayley = skew((R-I)*inv(I+R));

end

