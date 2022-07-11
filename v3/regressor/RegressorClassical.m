function [Y, W] = RegressorClassical( model, q, qd,qdd)
model = model.postProcessModel();
a_grav = model.getGravity();
if ~iscell(q)
    [q,qd,qdd] = confVecToCell(model,q,qd,qdd);

end

% Hotfix for symbolic REMOVE
Y = zeros(model.NV, model.N_RB*10);
%Y = sym(zeros(model.NV, model.N_RB*10)); % REMOVE
v = {};
a = {};
W = {};

for i = 1:model.NB
  vp = model.getParentVariable(i, v);
  ap = model.getParentVariable(i, a, -a_grav);
  [Xup{i}, S{i}, Sd{i}, v{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i}, qd{i}, vp);
  a{i} = Xup{i}*ap + Sd{i}*qd{i} + S{i}*qdd{i};
  F{i} = groupRegressor(a{i}, v{i});
end

for i = model.NB:-1:1
  ii = model.vinds{i};
  param_inds = model.param_inds{i};
  Y(ii,param_inds) = S{i}.' * F{i};
  W{i,i} = F{i};
  Fi = Xup{i}.' * F{i};
  
  j = model.parent(i);
  while j ~= 0
      jj = model.vinds{j};
      Y(jj,param_inds) = S{j}.' * Fi;
      W{j,i} = Fi;
      Fi = Xup{j}.' * Fi;
      j = model.parent(j);
  end
end

% pad zeros and concatenate W cells into matrix
for i = 2:size(W,1) 
    for j = 1:(i-1)
        W{i,j} = zeros( size(W{i,i},1), size(W{j,j},2) );
    end
end
W = cell2mat(W);

