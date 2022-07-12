%% Calculate Derivatives of the Fourier Trajectory w.r.t. params

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
%  dq           spline pos derivative vector                 
%  dqdot        (optional) spline vel derivative vector       
%  dqddot       (optional) spline acc derivative vector       

%% Implementation
function [dq, dqdot, dqddot] = getFourierDerivative(params, w, t)
    n = size(params,2);
    m = floor(size(params,1)/2); % k = 1,2,...m
    
    num_t = size(t,2);
    
    dq = zeros(n, 2*m, num_t);
    dqdot = zeros(n, 2*m, num_t);
    dqddot = zeros(n, 2*m, num_t);
    
    for i = 1:n
        for j = 1:num_t
            for k = 1:m
                dq(i,2*k-1,j)     =   sin(k*w*t(j));
                dq(i,2*k,j)       =   cos(k*w*t(j));
                dqdot(i,2*k-1,j)  =   k*w*cos(k*w*t(j));
                dqdot(i,2*k,j)    = - k*w*sin(k*w*t(j));
                dqddot(i,2*k-1,j) = - k*k*w*w*sin(k*w*t(j));                
                dqddot(i,2*k,j)   = - k*k*w*w*cos(k*w*t(j));                
            end
        end
    end
end