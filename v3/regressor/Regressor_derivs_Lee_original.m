%% Recursive Regressor Derivatives Calculator
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                        [Size]
%  A            i-th body screws from i-th frame     6*n    n = num of joints
%  M            initial relative frames M_{i,i-1}    4*4*n
%  q            joint angles                         n*1
%  V            spatial velocities                   6*n
%  Vdot         spatial accelerations                6*n
%  dq           joint angle derivative               n*m*n    m = num of parameters
%  dV           spatial velocitie derivative         6*n*m*n
%  dVdot        spatial acceleration derivative      6*n*m*n
%  W            (optional) sub-regressor matrix      6n*10n

%% Outputs
% [Name]  [Description]                                [Size]
%  dY      regressor derivative matrix                  n*10n*m*n
%  dW      (optional) sub-regressor derivative matrix   6n*10n*m*n                     

%% Implementation
function [dY, varargout] = getRegressorDerivativesRecursive(A,M,q,V,Vdot,dq,dV,dVdot,varargin)
    %% Initialization
    n      = size(q,1);         % number of joints
    m      = size(dq,2);        % number of parameters
    dW     = zeros(6*n, 10*n, m, n);
    
    if     nargin == 8          % no optional inputs
        [Y, W] = getRegressorRecursive(A,M,q,V,Vdot);
    elseif nargin == 9          % optional W
        W = varargin{1};
    else
        error('Init Error: undefined number of inputs');
    end
    
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
        
        for p=1:m
            for l = 1:n
                dW(6*i-5:6*i,10*i-9:10*i,p,l) = convertVelocityToRegressor(dVdot(:,i,p,l)) - small_ad(dV(:,i,p,l))'*convertVelocityToRegressor(V(:,i)) ...
                             - small_ad(V(:,i))'*convertVelocityToRegressor(dV(:,i,p,l));
                for k = i+1:n
                    dW(6*i-5:6*i,10*k-9:10*k,p,l) = Ad_T(:,:,i+1)'*dW(6*i+1:6*i+6,10*k-9:10*k,p,l) - Ad_T(:,:,i+1)'*small_ad(A(:,i+1))'*W(6*i+1:6*i+6,10*k-9:10*k)*dq(i+1,p,l);
                end                
            end
        end
    end
    
    dY = zeros(n, 10*n, m, n);
    for p =1:m
        for l = 1:n
            dY(:,:,p,l) = diagA'*dW(:,:,p,l);
        end
    end

    if nargout > 1
        varargout{1} = dW;
    end
end
