function  [tau,info] = ID_SlotineLI_rotoror( model, q, qd , qd_r, qdd, factorFunction)

if nargin == 5
    factorFunction = @(I,v)(factorFunctions(I,v));
end

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q)
    [q,qd,qd_r,qdd] = confVecToCell(model,q,qd,qd_r,qdd);
    
end

a_grav = get_gravity(model);
p = model.parent;
v = {};
w = {};
wd= {};

for i = 1:model.NB
  
  [Xup{i}, S{i}, Sd{i}, v{i}] = model.joint{i}.kinematics(q{i}, qd{i}, vp);
  wp = getParentVariable(model, i, w);
  wdp= getParentVariable(model, i, wd, -a_grav);
  
  v{i} = Xup{i}*vp + S{i}*qd{i};
  w{i} = Xup{i}*wp + S{i}*qd_r{i};
  wd{i}= Xup{i}*wdp + S{i} * qdd{i} + Sd{i}*qd_r{i};
  f{i} = model.I{i}*wd{i} + factorFunction(model.I{i}, v{i})*w{i};
end

tau = q{1}(1)*0 + zeros(model.NV,1);
for i = model.NB:-1:1
  ii = model.vinds{i};
  tau(ii) = S{i}' * f{i} ;
  if p(i) ~= 0
    f{p(i)} = f{p(i)} + Xup{i}'*f{i} ;
  end
end
