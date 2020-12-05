% There are three main function types created: 
% 1) Output functions: These provide some output, linear in the inertial
%    paramters. For instance, given some velocities v1, v2 one such output
%    gives y = v1'*I*v2 expressed linearly w.r.t. inertial paramters a.
%
% 2) Rate functions: These provide the rate of change in some parameters
%    as a function of an input velocity.
%
% 3) Transform functions: These propogate parameters across a link as a 
%    function of an input spatial transform.
%    
% Create a bunch of symbolic variables that will help us to write support
% functions using matlabFunction.
clear
syms m mcx mcy mcz Ixx Iyy Izz Iyz Ixz Ixy real
a  = [m mcx mcy mcz Ixx Iyy Izz Iyz Ixz Ixy]'; % Inertial parameter vector
a2 = sym('a2',[10 1],'real');

v1 = sym('v1',[6,1],'real');    vv1 = sym('vv1',[12,1],'real');      
v2 = sym('v2',[6,1],'real');    vv2 = sym('vv2',[12,1],'real');
     
X = sym('X',[6,6],'real');  % Mock spatial transform matrix
K = sym('Q',[6,6],'real');  % Mock matrix for outputs of form tr( K * I)
d = sym('d',[10 1],'real'); % Vector of dual inertial parameters

I = inertiaVecToMat(a);     % 6x6 spatial inertia
I2= inertiaVecToMat(a2); 

%% Create Functions for inertial parameters
y1= v1'*I*v2;          % Output for inner product between v1 and v2 
C1 = jacobian(y1,a);   % Linear relation: y1 = C1 * a
matlabFunction(C1,'File','Output_InnerProduct','vars', {v1,v2});

y2 = vv1'*[I zeros(6) ; zeros(6) I2]*vv2; % Likewise for motors
C2 = jacobian(y2,[a' a2']');              % Linear relation: y2 = C2 * [a ; a2] 
matlabFunction(C2,'File','Output_InnerProduct_rotor','vars', {vv1,vv2});

a_dot = inertiaMatToVec(crf(v1)*I - I*crm(v1)); % Rate of change in inertia paramters when moving with v1
A_phi = jacobian(a_dot , a);                    % Linear relation: a_dot = A_phi * a
matlabFunction(A_phi,'File','Rate_Parameters','vars',{v1});

a_trans = inertiaMatToVec(X' * I * X); % Parameters after transformation
A_X   = jacobian(a_trans, a);          % Linear relation: a_trans = A_X * a
matlabFunction(A_X,'File','Transform_Parameters','vars',{X});