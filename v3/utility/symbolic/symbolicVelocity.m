function [ param_vec, param_struct ] = symbolicVelocity( body_number )
% symbolicVelocity returns a symbolic spatial velocity vector 
%
% Input: (optional) body_number (let's call it i)
%         
%
% Ouput: param_vec = [omega_x_i omega_y_i omega_z_i v_x_i v_y_i v_z_i] (symbolic)
%
%        params_struct (same, but as structure with elements wx, wy,wz, vx,vy,vz etc.)

if nargin == 0
    append = '';
else
    append = ['_' num2str(body_number)];
end

wx  = sym(['omega_x'  append], 'real');
wy  = sym(['omega_y'  append], 'real');
wz  = sym(['omega_z'  append], 'real');
vx  = sym(['v_x'  append], 'real');
vy  = sym(['v_y'  append], 'real');
vz  = sym(['v_z'  append], 'real');


param_vec = [wx wy wz vx vy vz]';

param_struct.wx  = wx;
param_struct.wy  = wy;
param_struct.wz  = wz;
param_struct.vx  = vx;
param_struct.vy  = vy;
param_struct.vz  = vz;


end

