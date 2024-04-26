function [ param_vec, param_struct ] = symbolicForce( body_number )
% symbolicVelocity returns a symbolic spatial force vector 
%
% Input: (optional) body_number (let's call it i)
%         
%
% Ouput: param_vec = [n_x_i n_y_i n_z_i f_x_i f_y_i f_z_i] (symbolic)
%
%        params_struct (same, but as structure with elements nx, ny,nz, fx,fy,fz etc.)

if nargin == 0
    [n, n_struct] = symbolicCartesian('n');
    [f, f_struct] = symbolicCartesian('f');
else
    [n, n_struct] = symbolicCartesian('n',body_number);
    [f, f_struct] = symbolicCartesian('f',body_number);  
end

param_vec = [n ; f];
param_struct = n_struct;
entries = fields(f_struct);

for i = 1:length(entries)
    field_name = entries{i};
    param_struct.(field_name) = f_struct.(field_name);
end

end

