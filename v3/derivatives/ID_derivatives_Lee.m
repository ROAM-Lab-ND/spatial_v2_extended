function  [dtau, dV, dVd] = ID_derivatives_Lee( model, q, qd, qdd, dq, dqd, dqdd, F )

% derivatives of inverse dynamics w.r.t. trajectory parameters

% Adapted from Bryan Dongik Lee, 2018
% Original code at: https://github.com/SNURobotics/optimal-excitation/tree/master/libraries/dynamics

model = model.postProcessModel();
a_grav = model.getGravity();
m      = size(dq,2);        % number of parameters

if ~iscell(q)
    [q, qd, qdd] = confVecToCell(model,q,qd,qdd);
end

if ~iscell(dq)
    [dq, dqd, dqdd] = confdVecToCell(model,dq,dqd,dqdd);
end

% F should be output of ID, already as a cell array

dtau = {};
dV = {};
dVd = {};
v = {};
a = {};

for i = 1:model.NB
    vp = model.getParentVariable(i, v);
    ap = model.getParentVariable(i, a, -a_grav);
    [Xup{i}, S{i}, Sd{i}, v{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i}, qd{i}, vp);
    a{i} = Xup{i}*ap + Sd{i}*qd{i} + S{i}*qdd{i};

    if i==1 % for first joint
%         % note: these expressions are different than what is in the PhD thesis
%         dV{i} = S{i}*dqd{i};
%         dVd{i}= crm(dV{i})*S{i}*qd{i} + Sd{i}*dqd{i} + S{i}*dqdd{i};
        % expressions from PhD thesis:
        % TODO: could make more robust for if the first joint is an asbolute pair or triplet, not an issue now
        dV{i} = S{i}*dqd{i};
        dVd{i} = S{i}*dqdd{i} + Sd{i}*dqd{i} - crm(S{i})*dV{i}*qd{i};
    else
%         % note: these expressions are different than what is in the PhD thesis
%         dV{i} = Xup{i}*dV{i-1} - crm(S{i})*Xup{i}*v{i-1}*dq{i} + S{i}*dqd{i};
%         dVd{i} = Xup{i}*dVd{i-1} - crm(S{i})*Xup{i}*a{i-1}*dq{i} + crm(dV{i})*S{i}*qd{i} ...
%                     + Sd{i}*dqd{i} + S{i}*dqdd{i};
        % expressions from PhD thesis:
        dV{i} = Xup{i}*dV{i-1};
        dVd{i} = Xup{i}*dVd{i-1};
        for j=1:size(dq{i},1)
            dV{i} = dV{i} + S{i}(:,j)*dqd{i}(j,:) + crm(S{i}(:,j))*Xup{i}*v{i-1}*dq{i}(j,:);
        end
        for j=1:size(dq{i},1) % separate for loop since it depends on dV{i}
            dVd{i} = dVd{i} + S{i}(:,j)*dqdd{i}(j,:) + Sd{i}(:,j)*dqd{i}(j,:) - crm(S{i}(:,j))*Xup{i}*a{i-1}*dq{i}(j,:) - crm(S{i}(:,j))*dV{i}*qd{i}(j);
        end
%         dV{i} = S{i}*dqd{i} + Xup{i}*dV{i-1} + crm(S{i})*Xup{i}*v{i-1}*dq{i};
%         dVd{i} = S{i}*dqdd{i} + Xup{i}*dVd{i-1} + Sd{i}*dqd{i} - crm(S{i})*Xup{i}*a{i-1}*dq{i} - crm(S{i})*dV{i}*qd{i}; 
    end
end
  
dF = {};
dtau = {};

for i = model.NB:-1:1
    G{i} = model.I{i};
    % calculate dF
    if i==model.NB
        for p=1:m % for each parameter
            dF{i}(:,p) = G{i}*dVd{i}(:,p) + crf(dV{i}(:,p))*G{i}*v{i} + crf(v{i})*G{i}*dV{i}(:,p);
            % Original: dF(:,i,p,k) = G(:,:,i)*dVdot(:,i,p,k) - small_ad(dV(:,i,p,k))'*G(:,:,i)*V(:,i) - small_ad(V(:,i))'*G(:,:,i)*dV(:,i,p,k);  p=1:m,k=1:n
        end
    else
        for p=1:m
            dF{i}(:,p) = Xup{i+1}.'*dF{i+1}(:,p) + G{i}*dVd{i}(:,p) + crf(dV{i}(:,p))*G{i}*v{i} + crf(v{i})*G{i}*dV{i}(:,p);
            for j=1:size(dq{i+1},1)
                dF{i}(:,p) = dF{i}(:,p) - Xup{i+1}.'*crm(S{i+1}(:,j)).'*F{i+1}*dq{i+1}(j,p);
            end
            % Original: dF(:,i,p,k) = Ad_T(:,:,i+1)'* dF(:,i+1,p,k) - Ad_T(:,:,i+1)'*small_ad(A(:,i+1))'*F(:,i+1)*dq(i+1,p,k) + G(:,:,i)*dVdot(:,i,p,k) ...
            %                         - small_ad(dV(:,i,p,k))'*G(:,:,i)*V(:,i) - small_ad(V(:,i))'*G(:,:,i)*dV(:,i,p,k); 
        end
    end
    % calculate dtau
    dtau{i} = S{i}.'*dF{i};
    % Original: dtau(i,p,k) = A(:,i)'*dF(:,i,p,k);
end

end


%% function for converting trajectory derivatives to cell arrays of the correct sizes
function [dq_cell, dqd_cell, dqdd_cell] = confdVecToCell(model,dq,dqd,dqdd)
qj = 0;
vj = 0;

% TODO: sort out sizes!
dq_cell = cell(model.NB,1);
if nargin > 2
    dqd_cell = cell(model.NB,1);
end
if nargin > 3
    dqdd_cell = cell(model.NB,1);
end

for i = 1:model.NB
    
    dq_cell{i} = dq(qj+1:qj+model.joint{i}.nq,:);
    if nargin > 2
        dqd_cell{i} = dqd( vj+1 : vj+model.joint{i}.nv,:);
    end
    if nargin > 3
        dqdd_cell{i} = dqdd( vj+1 : vj+model.joint{i}.nv,:);
    end

    qj = qj+model.joint{i}.nq;
    vj = vj+model.joint{i}.nv;
end
end