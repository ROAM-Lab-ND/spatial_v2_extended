function Y = RegressorClassical( model, q, qd,qdd)
model = model.postProcessModel();
a_grav = model.getGravity();
if ~iscell(q)
    [q,qd,qdd] = confVecToCell(model,q,qd,qdd);

end

Y = zeros(model.NV, model.N_RB*10);
v = {};
a = {};

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
  Fi = Xup{i}.' * F{i};
  
  j = model.parent(i);
  while j ~= 0
      jj = model.vinds{j};
      Y(jj,param_inds) = S{j}.' * Fi;
      Fi = Xup{j}.' * Fi;
      j = model.parent(j);
  end
end
