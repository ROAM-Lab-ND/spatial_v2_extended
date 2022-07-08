%% Recursive Derivative Newton-Euler Inverse Dynamics Solver
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  A            i-th body screws from i-th frame   6*n    n = num of joints
%  M            initial relative frames M_{i,i-1}  4*4*n
%  q            joint angles                       n*1
%  qdot         joint vleocities                   n*1
%  G            link inertial matrices             6*6*n
%  V            spatial velocities                 6*n
%  Vdot         spatial accelerations              6*n
%  dq           dq/dp                              n*m*n,   m = num of parameters per joint
%  dqdot        dqdot/dp                           n*m*n
%  dqddot       dqddot/dp                          n*m*n
%  F            wrenches                           6*n
%  Vdot_0       (optional) base acceleration       6*1

%% Outputs
% [Name]  [Description]                                    [Size]
%  dtau    derivative of joint torques                      n*m*n
%  dV      (optional) derivative of spatial velocities      6*n*m*n
%  dVdot   (optional) derivative of spatial accelerations   6*n*m*n
%  dF      (optional) derivative of wrenches                6*n*m*n

%% Examples
% dtau = solveInverseDynamicsDerivatives(A,M,q,qdot,G,V,Vdot,dq,dqdot,dqddot,F)
% [dtau, dV, dVdot, dF] = solveInverseDynamicsDerivatives(A,M,q,qdot,G,V,Vdot,dq,dqdot,dqddot,F)

%% Implementation
function [dtau, varargout] = solveInverseDynamicsDerivatives(A,M,q,qdot,G,V,Vdot,dq,dqdot,dqddot,F,varargin)
    %% Initialization
    n      = size(q,1);         % number of joints
    m      = size(dq,2);        % number of parameters
    dV     = zeros(6,n,m,n);
    dVdot  = zeros(6,n,m,n);
    dF     = zeros(6,n,m,n);
    dtau   = zeros(n,m,n);
    
    V_0     = zeros(6,1);       % base velocity
    Vdot_0  = zeros(6,1);       % base acceleration
    dV_0    = zeros(6,1);       % derivative of base velocity
    dVdot_0 = zeros(6,1);       % derivative of base acceleration

    if     nargin == 11         % no optional inputs
    elseif nargin == 12         % optional base acceleration
        Vdot_0 = varargin{1};        
    else
        error('Init Error: undefined number of inputs');
    end
    
    T    = zeros(4,4,n); % T_{i,i-1}
    Ad_T = zeros(6,6,n); % Ad_T_{i,i-1}
    
    %% Forward Recursion
    for i = 1:n
        T(:,:,i)    = exp_se3(-A(:,i)*q(i))*M(:,:,i);
        Ad_T(:,:,i) = large_Ad(T(:,:,i));
        if i == 1
            for p = 1:m
                for k = 1:n
                    dV(:,i,p,k)    = Ad_T(:,:,i)*dV_0    - small_ad(A(:,i))*Ad_T(:,:,i)*V_0*dq(i,p,k)    + A(:,i)*dqdot(i,p,k);
                    dVdot(:,i,p,k) = Ad_T(:,:,i)*dVdot_0 - small_ad(A(:,i))*Ad_T(:,:,i)*Vdot_0*dq(i,p,k) + small_ad(dV(:,i,p))*A(:,i)*qdot(i) ...
                                    + small_ad(V(:,i))*A(:,i)*dqdot(i,p,k) + A(:,i)*dqddot(i,p,k);
                end
            end
        else
            for p = 1:m
                for k = 1:n
                    dV(:,i,p,k)    = Ad_T(:,:,i)*dV(:,i-1,p,k)    - small_ad(A(:,i))*Ad_T(:,:,i)*V(:,i-1)*dq(i,p,k)    + A(:,i)*dqdot(i,p,k);
                    dVdot(:,i,p,k) = Ad_T(:,:,i)*dVdot(:,i-1,p,k) - small_ad(A(:,i))*Ad_T(:,:,i)*Vdot(:,i-1)*dq(i,p,k) + small_ad(dV(:,i,p))*A(:,i)*qdot(i) ...
                                    + small_ad(V(:,i))*A(:,i)*dqdot(i,p,k) + A(:,i)*dqddot(i,p,k);
                end   
            end
        end
    end

%     %% Backward Recursion
%     for i = n:-1:1
%         if i == n
%             for p = 1:m
%                 for k = 1:n
%                     dF(:,i,p,k) = G(:,:,i)*dVdot(:,i,p,k) - small_ad(dV(:,i,p,k))'*G(:,:,i)*V(:,i) - small_ad(V(:,i))'*G(:,:,i)*dV(:,i,p,k);                
%                 end
%             end
%         else
%             for p = 1:m
%                 for k = 1:n
%                     dF(:,i,p,k) = Ad_T(:,:,i+1)'* dF(:,i+1,p,k) - Ad_T(:,:,i+1)'*small_ad(A(:,i+1))'*F(:,i+1)*dq(i+1,p,k) + G(:,:,i)*dVdot(:,i,p,k) ...
%                                 - small_ad(dV(:,i,p,k))'*G(:,:,i)*V(:,i) - small_ad(V(:,i))'*G(:,:,i)*dV(:,i,p,k);                
%                 end
%             end
%         end
%         for p = 1:m
%             for k = 1:n
%                 dtau(i,p,k) = A(:,i)'*dF(:,i,p,k);
%             end
%         end
%     end

    if nargout > 1
        varargout{1} = dV;
        varargout{2} = dVdot;
        varargout{3} = dF;
    end
end

