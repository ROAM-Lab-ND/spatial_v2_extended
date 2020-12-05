function [ quat ] = rotToQuat( R )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    angAx = rotToAngleAxis(R);
    quat = angleAxisToQuat(angAx);
end

