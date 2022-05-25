function  [tau,info] = ID_SlotineLi( model, q, qd , qd_r, qdd, factorFunction)

if nargin == 5
    factorFunction = @(I,v)(factorFunctions(I,v));
end

if ~isfield(model,'nq')
    model = model.postProcessModel();
end
if ~iscell(q)
    [q,qd,qd_r,qdd] = model.confVecToCell(q,qd,qd_r,qdd);
    
end

a_grav = model.getGravity();
p = model.parent;
v = {};
w = {};
wd= {};

for i = 1:model.NB
  vp = getParentVariable(model, i, v);
  wp = getParentVariable(model, i, w);
  wdp= getParentVariable(model, i, wd, -a_grav);
  
  [Xup{i}, S{i}, Sd{i}, v{i}] = model.joint{i}.kinematics(q{i}, qd{i}, vp);
  
  w{i} = Xup{i}*wp + S{i}*qd_r{i};
  wd{i}= Xup{i}*wdp + S{i} * qdd{i} + Sd{i}*qd_r{i};
  f{i} = model.I{i}*wd{i} + factorFunction(model.I{i}, v{i})*w{i};
end

tau = q{1}(1)*0 + zeros(model.NV,1);
for i = model.NB:-1:1
  ii = model.vinds{i};
  tau(ii) = S{i}.' * f{i} ;
  if p(i) ~= 0
    f{p(i)} = f{p(i)} + Xup{i}.'*f{i} ;
  end
end
