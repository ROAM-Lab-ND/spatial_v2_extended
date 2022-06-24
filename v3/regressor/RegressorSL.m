function  [Y] = RegressorSL( model, q, qd,qd_r,qdd, factorFunction)

if nargin == 5
    factorFunction = @(I,v)(factorFunctions(I,v));
end
model = model.postProcessModel();

if ~iscell(q)
    [q,qd,qd_r,qdd] = confVecToCell(model,q,qd,qd_r,qdd);
    
end

Y = zeros(model.NV, model.N_RB*10);
a_grav = model.getGravity();

v = {};
w = {};
wd= {};

for i = 1:model.NB
  vp = getParentVariable(model, i, v);
  wp = getParentVariable(model, i, w);
  wdp= getParentVariable(model, i, wd, -a_grav);
  
  [Xup{i}, S{i}, Sd{i}, v{i}] = model.joint{i}.kinematics(model.Xtree{i},q{i}, qd{i}, vp);
  
  w{i} = Xup{i}*wp + S{i}*qd_r{i};
  wd{i}= Xup{i}*wdp + S{i} * qdd{i} + Sd{i}*qd_r{i};
  
  F{i} = groupRegressor(wd{i}, v{i}, w{i}, factorFunction);
end

for i = model.NB:-1:1
  ii = model.vinds{i};
  param_inds = model.param_inds{i};
  Y(ii,param_inds) = S{i}' * F{i};
  Fi = Xup{i}' * F{i};
 
  j = model.parent(i);
  while j ~= 0
      jj = model.vinds{j};
      Y(jj,param_inds) = S{j}' * Fi;
      Fi = Xup{j}' * Fi;
      j = model.parent(j);
  end
end