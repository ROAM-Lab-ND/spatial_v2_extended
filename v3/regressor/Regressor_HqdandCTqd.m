function  [Y_Hqd, Y_CTqd]  = Regressor_HqdandCTqd( model, q , qd)
% Indirect regressors term  YHqd a = H qd and YCTqd = C^T qd

model = model.postProcessModel();
if ~iscell(q)
    [q, qd] = confVecToCell(model,q,qd);
end

Y_Hqd = zeros(model.NV, 10*model.N_RB);
Y_CTqd = zeros(model.NV, 10*model.N_RB);

v = {};
for i = 1:model.NB
  vp = model.getParentVariable(i, v);
  [Xup{i}, S{i}, Sd{i}, v{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i}, qd{i}, vp);
  
  hi = groupRegressor(v{i}, 0*v{i});
  param_inds = model.param_inds{i};
  
  j = i;
  while j > 0
    jj = model.vinds{j};
    Y_Hqd(jj, param_inds) = S{j}.' * hi;
    Y_CTqd(jj, param_inds)= Sd{j}.'* hi;
    hi = Xup{j}.' * hi;
	j = model.parent(j);
  end
end