%% Fourier Series Trajectory Generator a_k*sin(kwt) + b_k*cos(kwt)

% Adapted from Bryan Dongik Lee, 2018
% Original code at: https://github.com/SNURobotics/optimal-excitation/tree/master/libraries/trajectory

% helper for testing arm regressor and regressor derivative functions

%% Inputs
% [Name]       [Description]
%  params       spline coefficients               
%  w            base frequency (sin wt)
%  t            times to sample 

%% Outputs
% [Name]       [Description]
%  q            spline pos vector
%  qdot         (optional) spline vel vector
%  qddot        (optional) spline acc vector

%% Implementation
function [q, qdot, qddot] = makeFourier(params, w, t, varargin)
    n = size(params,2);
    m = floor(size(params,1)/2); % k = 1,2,...m
    
    num_t = size(t,2);
    
    q = zeros(n, num_t);        
    qdot = zeros(n, num_t);
    qddot = zeros(n, num_t);
    
    if nargin > 3
        q = q + varargin{1};
    end
    
    for i = 1:n
        for j = 1:num_t
            for k = 1:m
                q(i,j) = q(i,j) + params(2*k-1,i)*sin(k*w*t(j)) + params(2*k,i)*cos(k*w*t(j));
                qdot(i,j) = qdot(i,j) + k*w*params(2*k-1,i)*cos(k*w*t(j)) - k*w*params(2*k,i)*sin(k*w*t(j));
                qddot(i,j) = qddot(i,j) - k*k*w*w*params(2*k-1,i)*sin(k*w*t(j)) - k*k*w*w*params(2*k,i)*cos(k*w*t(j));                
            end
        end
    end
end