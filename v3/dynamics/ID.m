function  [tau, out] = ID( model, q, qd, qdd, f_ext )

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

a_grav = get_gravity(model);
if ~iscell(q)
    [q, qd, qdd] = confVecToCell(model,q,qd,qdd);
end


for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    a{i} = Xup{i}*(-a_grav) + S{i}*qdd{i};
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd{i} + crm(v{i})*vJ;
  end
  h{i} = model.I{i}*v{i};
  f{i} = model.I{i}*a{i} + crf(v{i})*h{i};
  
  % Extra data for rotors
  if model.has_rotor(i)
      [ XJ_rotor, S_rotor{i} ] = jcalc( model.jtype_rotor{i}, q{i}*model.gr{i} );
      S_rotor{i} = S_rotor{i} * model.gr{i};
      vJ_rotor = S_rotor{i} * qd{i};
      Xup_rotor{i} = XJ_rotor * model.Xrotor{i};
      if model.parent(i) == 0
          v_rotor{i} = vJ_rotor;
          a_rotor{i} = Xup_rotor{i}*(-a_grav) + S_rotor{i}*qdd{i};
      else
          v_rotor{i} = Xup_rotor{i}*v{model.parent(i)} + vJ_rotor;
          a_rotor{i} = Xup_rotor{i}*a{model.parent(i)} + S_rotor{i}*qdd{i} + crm(v_rotor{i})*vJ_rotor;
      end
      f_rotor{i} = model.I_rotor{i}*a_rotor{i} + crf(v_rotor{i})*model.I_rotor{i}*v_rotor{i};
  end
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
      
  if model.has_rotor(i) % Modified backward pass
      tau(ii) = tau(ii) + S_rotor{i}' * f_rotor{i};
      if p ~= 0
        f{p} = f{p} + Xup_rotor{i}.'*f_rotor{i};
      end
  end
end

out.Xup = Xup;
out.v = v;
out.h = h;
out.a = a;
out.f = f;
out.tau = tau;