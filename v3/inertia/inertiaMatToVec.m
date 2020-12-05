function [ a ] = inertiaMatToVec( I )
% inertiaMatToVec converts a 6x6 spatial inertia matrix into a vector of 
% 10 inertial parameters 
%
% Input: I (6x6) spatial inertia matrix
%
% Outpt: a (10x1) Vector of inertial parameters
%         a = [m hx hy hz Ixx Iyy Izz Ixy Ixz Iyz]'
%

    h = skew(I(1:3,4:6));
    a = [I(6,6)  h(1) h(2) h(3) I(1,1) I(2,2) I(3,3) I(3,2) I(3,1) I(2,1) ]';
end

