function  qdd = FDab( model, q, qd, tau, f_ext )

% FDab  Forward Dynamics via Articulated-Body Algorithm
% FDab(model,q,qd,tau,f_ext,grav_accn)  calculates the forward dynamics of
% a kinematic tree via the articulated-body algorithm.  q, qd and tau are
% vectors of joint position, velocity and force variables; and the return
% value is a vector of joint acceleration variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.  The format of f_ext is
% explained in the source code of apply_external_forces.


if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q) || ~iscell(qd)
    [q, qd, tau] = confVecToCell(model,q,qd,tau);
end

symzero = q{1}(1)*0;
a_grav = model.getGravity();

v = {};
a = {};

for i = 1:model.NB
  vp = model.getParentVariable(i, v); 
  [Xup{i}, S{i}, Sd{i}, v{i}] = model.joint{i}.kinematics(model.Xtree{i},q{i}, qd{i}, vp);
  c{i} = Sd{i}*qd{i};
  IA{i} = symzero+model.I{i};
  pA{i} = crf(v{i}) * model.I{i} * v{i};
end

if nargin == 5
  pA = apply_external_forces( model.parent, Xup, pA, f_ext );
end

for i = model.NB:-1:1
  p = model.parent(i);
  
  U{i}  = IA{i} * S{i};
  d{i}  = S{i}.' * U{i};

  u{i} = tau{i} - S{i}.'*pA{i} - U{i}.'*c{i};
  U{i} = Xup{i}.' * U{i};

  if model.parent(i) ~= 0
    Ia = Xup{i}.' * IA{i} * Xup{i} - U{i}/d{i}*U{i}.';
    pa = Xup{i}.' *(pA{i} + IA{i}*c{i}) + U{i} * (d{i}\u{i});

    IA{p} = IA{p} + Ia;
    pA{p} = pA{p} + pa;
  end
end

qdd = q{1}(1)*0 + zeros(model.NV,1);
for i = 1:model.NB
  ap = getParentVariable(model, i, a, -a_grav);
  ii = model.vinds{i};
  qdd(ii) = d{i}\(u{i} - U{i}.'*ap);
  a{i} = Xup{i} * ap + S{i}*qdd(ii) + c{i};
end