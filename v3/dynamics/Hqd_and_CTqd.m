function  [Hqd, CTqd]  = Hqd_and_CTqd( model, q , qd)

if ~isfield(model,'nq')
    model = model.postProcessModel();
end
if ~iscell(q)
    [q, qd] = model.confVecToCell(q,qd);
end

v = {};
for i = 1:model.NB
  vp = model.getParentVariable(i, v);
  [Xup{i}, S{i}, Sd{i}, v{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i}, qd{i}, vp); 
  h{i} = model.I{i}*v{i};
end

Hqd = q{1}(1)*0 + zeros(model.NV,1);
CTqd= q{1}(1)*0 + zeros(model.NV,1);

for i = model.NB:-1:1   
  ii = model.vinds{i};
  Hqd(ii) = S{i}.' * h{i};
  CTqd(ii)= Sd{i}.'* h{i};
  if model.parent(i) ~= 0
     p = model.parent(i);
     h{p} = h{p} + Xup{i}.'*h{i};
  end
end

