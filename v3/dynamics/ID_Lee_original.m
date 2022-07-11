%% Recursive Newton-Euler Inverse Dynamics Solver
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  A            i-th body screws from i-th frame   6*n
%  M            initial relative frames M_{i,i-1}  4*4*n
%  q            joint angles                       n*1
%  qdot         joint vleocities                   n*1
%  qddot        joint accelerations                n*1
%  G            link inertial matrices             6*6*n
%  (optional1) and (optional2)
%  Vdot_0       (optional1) base acceleration      6*1
%  F_desired    (optional2) end-effector force     6*1
%  T_end        (optional2) end-effector frame     4*4

%% Outputs
% [Name]  [Description]                      [Size]
%  tau     joint torques                      n*1
%  V       (optional) spatial velocities      6*n
%  Vdot    (optional) spatial accelerations   6*n
%  F       (optional) wrenches                6*n

%% Examples
% tau = solveInverseDynamics(A,M,q,qdot,qddot,G)
% [tau, V, Vdot] = solveInverseDynamics(A,M,q,qdot,qddot,G,Vdot_0,F_desired,T_end)
% [tau, V, Vdot, F] = solveInverseDynamics(A,M,q,qdot,qddot,G,Vdot_0)

%% Implementation
function [tau, varargout] = ID_Lee_original(A,M,q,qdot,qddot,G,varargin)
    %% Initialization
    n     = size(q,1);          % number of joints
    V     = zeros(6,n);
    Vdot  = zeros(6,n);
    F     = zeros(6,n);
    tau   = zeros(n,1);
    
    V_0    = zeros(6,1);        % base velocity
    Vdot_0 = zeros(6,1);        % base acceleration
    F_tip  = zeros(6,1);        % end-effector force
    T_end  = eye(4,4);          % end-effector pose T_{i+1,i}

    if     nargin == 6          % no optional inputs
    elseif nargin == 7          % optional input 1
        Vdot_0 = varargin{1};        
    elseif nargin == 8          % optional input 2
        F_tip  = varargin{1};
        T_end  = varargin{2};
    elseif nargin == 9          % optional input 1 & 2
        Vdot_0 = varargin{1};    
        F_tip  = varargin{2};
        T_end  = varargin{3};
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
            V(:,i)    = Ad_T(:,:,i)*V_0    + A(:,i)*qdot(i);
            Vdot(:,i) = Ad_T(:,:,i)*Vdot_0 + small_ad(V(:,i))*A(:,i)*qdot(i) + A(:,i)*qddot(i);
        else
            V(:,i)    = Ad_T(:,:,i)*V(:,i-1)    + A(:,i)*qdot(i);
            Vdot(:,i) = Ad_T(:,:,i)*Vdot(:,i-1) + small_ad(V(:,i))*A(:,i)*qdot(i) + A(:,i)*qddot(i);
        end
    end

    %% Backward Recursion
    for i = n:-1:1
        if i == n
            F(:,i) = large_Ad(T_end)'*F_tip + G(:,:,i)*Vdot(:,i) - small_ad(V(:,i))'*G(:,:,i)*V(:,i);
        else
            F(:,i) = Ad_T(:,:,i+1)'*F(:,i+1) + G(:,:,i)*Vdot(:,i) - small_ad(V(:,i))'*G(:,:,i)*V(:,i);
        end
        tau(i) = A(:,i)'*F(:,i);
    end
    
    if nargout > 1
        varargout{1} = V;
        varargout{2} = Vdot;
        varargout{3} = F;    
    end
end

