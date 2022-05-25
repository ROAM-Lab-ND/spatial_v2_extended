function  [q_dot, p_dot]  = HamiltonianDynamics( model, q , p, tau, f_ext)

if ~isfield(model,'nq')
    model = model.postProcessModel();
end

model_no_grav = model;
model_no_grav.gravity = [0 0 0]';

qd    = FDab( model_no_grav , q, zeros(model.NV,1), p);

if nargin == 5
    tau_g = ID(model, q, zeros(model.NV,1), zeros(model.NV,1), f_ext);
else
    tau_g = ID(model, q, zeros(model.NV,1), zeros(model.NV,1));
end

if ~iscell(q)
    [q, qd, p] = model.confVecToCell(q,qd,p);
end

v = {};
for i = 1:model.NB
  vp = model.getParentVariable(i, v);
  [Xup{i}, S{i}, Sd{i}, v{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i}, qd{i}, vp); 
  h{i} = model.I{i}*v{i};
end

CTqd= q{1}(1)*0 + zeros(model.NV,1);
q_dot = q{1}(1)*0 + zeros(model.NV,1);

for i = model.NB:-1:1   
  ii = model.vinds{i};
  q_dot(ii) = qd{i};
  CTqd(ii)= Sd{i}.'* h{i};
  if model.parent(i) ~= 0
     p = model.parent(i);
     h{p} = h{p} + Xup{i}.'*h{i};
  end
end

p_dot = CTqd + tau - tau_g;