function [ cayley ] = rpyToCayley( rpy )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    cayley = rotToCayley(rpyToRot(rpy));
end

