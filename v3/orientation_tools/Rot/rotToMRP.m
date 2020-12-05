function [ mrp ] = rotToMRP( R )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    angAx = rotToAngleAxis(R);
    [ang ax] = angleAxisToAngleAndAxis(angAx);
    mrp = [tan(ang/4) * ax];
end

