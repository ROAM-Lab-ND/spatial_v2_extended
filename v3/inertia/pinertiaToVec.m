function [ a ] = pinertiaToVec( J )
% pinertiaToVec converts a 4x4 Pseudo-Inertia matrix to a vector of 10
% inertial parameters. 
%
% Input:  J (4x4) psuedo inertia
% Output: a (10x1) Vector of inertial parameters
%         a = [m hx hy hz Ixx Iyy Izz Ixy Ixz Iyz]'
%
% The pseudo-inertia is discussed further in
%   Linear Matrix Inequalities for Physically Consistent Inertial Parameter 
%   Identification: A Statistical Perspective on the Mass Distribution, by 
%   Wensing, Kim, Slotine

    m = J(4,4);
    h = J(1:3,4);
    c = h/m;
    E = J(1:3,1:3);
    Ibar = trace(E)*eye(3) - E;
    
    a = inertiaMatToVec([Ibar skew(h) ; skew(h)' m*eye(3)]);
end

