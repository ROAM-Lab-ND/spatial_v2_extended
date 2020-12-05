function [ Pinertia ] = inertiaMatToPinertia( I )
% inertiaMatToPinertia converts a 6x6 spatial inertia matrix to a 4x4
% pseudo-inertia matrix
%
% Input:  I (6x6) spatial inertia matrix
% Output: J (4x4) Pseudo-inertia matrix
% 
% The pseudo-inertia is discussed further in
%   Linear Matrix Inequalities for Physically Consistent Inertial Parameter 
%   Identification: A Statistical Perspective on the Mass Distribution, by 
%   Wensing, Kim, Slotine

    h = skew(I(1:3,4:6));
    Ibar = I(1:3,1:3);
    m = I(6,6);
    Pinertia = [ 1/2*trace(Ibar)*eye(3)-Ibar h ; h' m ];
end

