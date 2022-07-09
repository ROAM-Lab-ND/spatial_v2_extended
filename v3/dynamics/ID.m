function  [tau, out] = ID( model, q, qd, qdd, f_ext )

% ID  Inverse Dynamics via Recursive Newton-Euler Algorithm
% ID(model,q,qd,qdd,f_ext) calculates the inverse dynamics of a kinematic
% tree via the recursive Newton-Euler algorithm.  q, qd and qdd are vectors
% of joint position, velocity and acceleration variables; and the return
% value is a vector of joint force variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.  The format of f_ext is
% explained in the source code of apply_external_forces.


model = model.postProcessModel();
a_grav = model.getGravity();

if ~iscell(q)
    [q, qd, qdd] = confVecToCell(model,q,qd,qdd);
end




v = {};
a = {};

for i = 1:model.NB
  vp = model.getParentVariable(i, v);
  ap = model.getParentVariable(i, a, -a_grav);
  
  [Xup{i}, S{i}, Sd{i}, v{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i}, qd{i}, vp);

  a{i} = Xup{i}*ap + Sd{i}*qd{i} + S{i}*qdd{i};
  h{i} = model.I{i}*v{i};
  f{i} = model.I{i}*a{i} + crf(v{i})*h{i};
end

out.f0 = f;

if nargin == 5
  f = apply_external_forces( model.parent, Xup, f, f_ext );
end

% This line ensures that ID is symbolic or Casadi compatible
tau = q{1}(1)*0 + zeros(model.NV,1);
    
for i = model.NB:-1:1
  p = model.parent(i);  
  ii = model.vinds{i};
  tau(ii) = S{i}.' * f{i};
  if p ~= 0
        f{p} = f{p} + Xup{i}.'*f{i} ;
  end
end

out.Xup = Xup;
out.S = S;
out.Sd = Sd;
out.v = v;
out.h = h;
out.a = a;
out.f = f;
out.tau = tau;