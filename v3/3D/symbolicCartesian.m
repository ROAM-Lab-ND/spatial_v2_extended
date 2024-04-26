function [ param_vec, param_struct ] = symbolicCartesian( symbol, body_number )
% symbolicCartesian returns a symbolic cartesian vector 
%
% Input: (optional) body_number (let's call it i)
%         
%
% Ouput: param_vec = [symbol_x_i symbol_y_i symbol_z_i] (symbolic)
%
%        params_struct (same, but as structure with elements symbolx, symboly,symbolz)

if nargin <= 1
    append = '';
else
    append = ['_' num2str(body_number)];
end

sx  = sym([symbol '_x'  append], 'real');
sy  = sym([symbol '_y'  append], 'real');
sz  = sym([symbol '_z'  append], 'real');

param_vec = [sx sy sz]';

param_struct.([symbol 'x'])  = sx;
param_struct.([symbol 'y'])  = sy;
param_struct.([symbol 'z'])  = sz;

end

