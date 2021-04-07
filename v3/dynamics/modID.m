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

if sum(model.has_rotor) > 1
    error('modID does not support rotors');
end

a_grav = get_gravity(model);
if ~iscell(q)
    [q, qd, qdd, lambda] = confVecToCell(model,q,qd,qdd, lambda);
end

out=q{1}(1)*0 + 0;
for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  wJ = S{i}*lambda{i};
  
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    w{i} = wJ;
    a{i} = Xup{i}*(-a_grav) + S{i}*qdd{i};
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    w{i} = Xup{i}*w{model.parent(i)} + wJ;
    a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd{i} + crm(v{i})*vJ;
  end
  out = out + w{i}.'*(model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i});
  
end