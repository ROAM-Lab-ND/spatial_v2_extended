function [ I ] = inertiaVecToMat( a )
% inertiaVecToMat converts a vector of 10 inertial parameters to 6x6 spatial
% inertia matrix
%
% Input: a (10x1) Vector of inertial parameters
%         a = [m hx hy hz Ixx Iyy Izz Ixy Ixz Iyz]'
%
% Ouput: I (6x6) spatial inertia matrix
%

Ibar = [a(5)  a(10) a(9) ; 
        a(10) a(6)  a(8) ; 
        a(9)  a(8)  a(7) ];
h    = [a(2)  a(3)  a(4)]';
m    = a(1);
I = [Ibar skew(h) ; skew(h)' m*eye(3)];

end

