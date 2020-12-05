function [ Pinertia ] = inertiaVecToPinertia( a )
% pinertiaToVec converts a vector of 10 inertial parameters to 4x4 
% Pseudo-Inertia matrix.
%
% Input:  J (4x4) psuedo inertia
% Output: a (10x1) Vector of inertial parameters
%         a = [m hx hy hz Ixx Iyy Izz Ixy Ixz Iyz]'
%
% The pseudo-inertia is discussed further in:
%  Linear Matrix Inequalities for Physically Consistent Inertial Parameter 
%  Identification: A Statistical Perspective on the Mass Distribution, by 
%  Wensing, Kim, Slotine

    I = inertiaVecToMat(a);
    h = skew(I(1:3,4:6));
    Ibar = I(1:3,1:3);
    m = I(6,6);
    Pinertia = [ 1/2*trace(Ibar)*eye(3)-Ibar h ; h' m ];
end

