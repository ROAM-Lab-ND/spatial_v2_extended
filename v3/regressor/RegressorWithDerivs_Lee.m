function [Y, dY, W, dW] = RegressorWithDerivs_Lee( model, q, qd, qdd, dq, dqd, dqdd)

model = model.postProcessModel();
m      = size(dq,2);        % number of parameters
% a_grav = model.getGravity();

if ~iscell(q)
    [q,qd,qdd] = confVecToCell(model,q,qd,qdd);
end

if ~iscell(dq)
    [dq, dqd, dqdd] = confdVecToCell(model,dq,dqd,dqdd);
end

Wcell = {};
dWcell = {};
W = [];
dW = [];

v = {};
a = {};

% call inverse dynamics and pull out necessary quantities
[tau,out] = ID(model,q,qd,qdd);
Xup = out.Xup;
v = out.v;
a = out.a;
S = out.S;

% call derivs of inverse dynamics
[dtau, dV, dVd] = ID_derivatives_Lee( model, q, qd, qdd, dq, dqd, dqdd, out.f );

% for i = 1:model.NB
%   vp = model.getParentVariable(i, v);
%   ap = model.getParentVariable(i, a, -a_grav);
%   [Xup{i}, S{i}, Sd{i}, v{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i}, qd{i}, vp);
%   a{i} = Xup{i}*ap + Sd{i}*qd{i} + S{i}*qdd{i};
% end

% asemble S_T block diagonal matrix
diagS = [];
for i=1:model.NB
    diagS = blkdiag(diagS, S{i});
end 

% recursion to calculate W
for i = model.NB:-1:1
    Wcell{i} = zeros(size(a{i},1),model.N_RB*10);
    dWcell{i} = zeros(size(a{i},1),model.N_RB*10, m);
    Wcell{i}(:,model.param_inds{i}) = vel_reg_Lee(a{i}) - crm(v{i}).'*vel_reg_Lee(v{i});
    for p = 1:m
        dWcell{i}(:,model.param_inds{i},p) = vel_reg_Lee(dVd{i}(:,p)) + crf(dV{i}(:,p))*vel_reg_Lee(v{i}) ...
                                                    + crf(v{i})*vel_reg_Lee(dV{i}(:,p));
    end
    for k = i+1:model.NB
        Wcell{i}(:,model.param_inds{k}) = Xup{i}.'*Wcell{i+1}(:,model.param_inds{k});
        for p = 1:m
            dWcell{i}(:,model.param_inds{k},p) = Xup{i}.'*dWcell{i+1}(:,model.param_inds{k},p);
            for j = 1:size(S{k},2)
                % TODO: need to check indices again here
                dWcell{i}(:,model.param_inds{k},p) = dWcell{i}(:,model.param_inds{k},p) - Xup{i}.'*crm(S{i+1}(:,j)).'*Wcell{i+1}(:,model.param_inds{k})*dq{k}(j,p);
            end
        end
    end
    W = [Wcell{i};W];
    dW = [dWcell{i};dW];
end

% convert from W to Y
Y = diagS.'*W;
dY = zeros(size(Y,1),size(Y,2),m);
for p=1:m % for each parameter
    dY(:,:,p) = diagS.'*dW(:,:,p);
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
