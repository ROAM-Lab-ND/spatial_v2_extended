function  [H,C,info] = HandC( model, q, qd, f_ext )

% HandC  Calculate coefficients of equation of motion
% [H,C]=HandC(model,q,qd,f_ext)  calculates the coefficients of the
% joint-space equation of motion, tau=H(q)qdd+C(d,qd,f_ext), where q, qd
% and qdd are the joint position, velocity and acceleration vectors, H is
% the joint-space inertia matrix, C is the vector of gravity,
% external-force and velocity-product terms, and tau is the joint force
% vector.  Algorithm: recursive Newton-Euler for C, and
% Composite-Rigid-Body for H.  f_ext is an optional argument specifying the
% external forces acting on the bodies.  It can be omitted if there are no
% external forces.  The format of f_ext is explained in the source code of
% apply_external_forces.

a_grav = model.getGravity();

if ~isfield(model,'nq')
    model = model.postProcessModel();
end
if ~iscell(q) || ~iscell(qd)
    [q, qd] = confVecToCell(model,q,qd);
end

v = {};
avp = {};

for i = 1:model.NB
  vp    = getParentVariable(model, i, v);
  avp_p = getParentVariable(model, i, avp, -a_grav);
  
  [Xup{i}, S{i}, Sd{i}, v{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i}, qd{i}, vp);
  avp{i} = Xup{i}*avp_p + Sd{i}*qd{i};
  fvp{i} = model.I{i}*avp{i} + crf(v{i})*model.I{i}*v{i};
end

if nargin == 4
  fvp = apply_external_forces( model.parent, Xup, fvp, f_ext );
end

C = q{1}(1)*0 + zeros(model.NV,1);
H = q{1}(1)*0 + zeros(model.NV);
for i =1:model.NB
    IC{i} = q{1}(1)*0 + model.I{i};
end

p = model.parent;
for i = model.NB:-1:1
  fh = IC{i} * S{i};
  ii = model.vinds{i};
  H(ii,ii) = S{i}.' * fh;
  C(ii,1) = S{i}.' * fvp{i};
  
  fh = Xup{i}.'* fh;
   
  j = i;
  while p(j) > 0
    j = p(j);
    jj = model.vinds{j};
    H(jj,ii) = S{j}.' * fh;
    H(ii,jj) = H(jj,ii).';
    fh = Xup{j}.' * fh;
  end
  
  if p(i) ~= 0
    IC{p(i)}  = IC{p(i)} + Xup{i}.'*IC{i}*Xup{i};
    fvp{p(i)} = fvp{p(i)} + Xup{i}.'*fvp{i};
  end
end

info.IC = IC;
info.fvp = fvp;
info.avp = avp;
info.Xup = Xup;
info.v = v;