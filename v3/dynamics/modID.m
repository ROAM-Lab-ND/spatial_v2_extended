function  out = modID( model, q, qd, qdd, lambda )

% ID  Inverse Dynamics via Recursive Newton-Euler Algorithm
% ID(model,q,qd,qdd,f_ext) calculates the inverse dynamics of a kinematic
% tree via the recursive Newton-Euler algorithm.  q, qd and qdd are vectors
% of joint position, velocity and acceleration variables; and the return
% value is a vector of joint force variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.  The format of f_ext is
% explained in the source code of apply_external_forces.

if ~isfield(model,'nq')
    model = postProcessModel(model);
end

a_grav = model.getGravity(model);
if ~iscell(q)
    [q, qd, qdd, lambda] = model.confVecToCell(q,qd,qdd, lambda);
end

out=q{1}(1)*0 + 0;
w = {};
v = {};
a = {};

for i = 1:model.NB
  vp = model.getParentVariable(i, v);
  ap = model.getParentVariable(i, a, -a_grav);
  wp = model.getParentVariable(i, w);
  [Xup{i}, S{i}, Sd{i}, v{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i}, qd{i}, vp);
  w{i} = Xup{i}*wp + S{i}*lambda{i};
  a{i} = Xup{i}*ap + S{i}*qdd{i} + Sd{i}*qd{i};
  out = out + w{i}.'*(model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i});
end