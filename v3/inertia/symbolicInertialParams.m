function [ param_vec, param_struct ] = symbolicInertialParams( body_number )
% symbolicInertialParams returns a symbolic vector of 10 inertial parameters 
%
% Input: (optional) body_number (let's call it i)
%         
%
% Ouput: param_vec = [m_i h_x_i h_y_i h_z_i I_xx_i I_yy_i I_zz_i I_xy_i I_xz_i
% I_yz_i]' (symbolic)
%
%        params_struct (same, but as structure with elements m, hx, etc.)

if nargin == 0
    append = '';
else
    append = ['_' num2str(body_number)];
end

m   = sym(['m'    append], 'real');
hx  = sym(['h_x'  append], 'real');
hy  = sym(['h_y'  append], 'real');
hz  = sym(['h_z'  append], 'real');
Ixx = sym(['I_xx' append], 'real');
Iyy = sym(['I_yy' append], 'real');
Izz = sym(['I_zz' append], 'real');
Iyz = sym(['I_yz' append], 'real');
Ixz = sym(['I_xz' append], 'real');
Ixy = sym(['I_xy' append], 'real');

param_vec = [m hx hy hz Ixx Iyy Izz Iyz Ixz Ixy]';

param_struct.m   = m;
param_struct.hx  = hx;
param_struct.hy  = hy;
param_struct.hz  = hz;
param_struct.Ixx = Ixx;
param_struct.Iyy = Iyy;
param_struct.Izz = Izz;
param_struct.Iyz = Iyz;
param_struct.Ixz = Ixz;
param_struct.Ixy = Ixy;

end

