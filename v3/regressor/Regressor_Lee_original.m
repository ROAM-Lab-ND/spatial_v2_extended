
%% Recursive Regressor Calculator
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  A            i-th body screws from i-th frame   6*n    n = num of joints
%  M            initial relative frames M_{i,i-1}  4*4*n
%  q            joint angles                       n*1
%  V            spatial velocities                 6*n
%  Vdot         spatial accelerations              6*n

%% Outputs
% [Name]  [Description]                                     [Size]
%  Y       regressor matrix (tau = Y*PHI)                    n*10n
%  W       (optional) sub-regressor matrix (Y = diag(A)'*W)  6n*10n                       

%% Implementation
function [Y, varargout] = getRegressorRecursive(A,M,q,V,Vdot)
    %% Initialization
    n      = size(q,1);         % number of joints  
    W      = zeros(6*n, 10*n);
    diagA  = zeros(6*n,n);
    for i=1:n
        diagA(6*i-5:6*i,i) = A(:,i);
    end  
    
    T    = zeros(4,4,n); % T_{i,i-1}
    Ad_T = zeros(6,6,n); % Ad_T_{i,i-1}
    
    %% Recursion
    for i = n:-1:1
        T(:,:,i)    = exp_se3(-A(:,i)*q(i))*M(:,:,i);
        Ad_T(:,:,i) = large_Ad(T(:,:,i));
        
        W(6*i-5:6*i,10*i-9:10*i) = convertVelocityToRegressor(Vdot(:,i)) - small_ad(V(:,i))'*convertVelocityToRegressor(V(:,i));
        for k = i+1:n
            W(6*i-5:6*i,10*k-9:10*k) = Ad_T(:,:,i+1)'*W(6*i+1:6*i+6,10*k-9:10*k);
        end
    end
    
    Y = diagA'*W;
    
    if nargout > 1
        varargout{1} = W;
    end
end

