function  qdd = FDab_rotor( model, q, qd, tau, f_ext )

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
a_grav = get_gravity(model);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    c{i} = zeros(size(a_grav));		% spatial or planar zero vector
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    c{i} = crm(v{i}) * vJ;
  end
  IA{i} = model.I{i};
  pA{i} = crf(v{i}) * model.I{i} * v{i};
  
  % Special stuff for rotors
  if model.has_rotor(i)
      [ XJ_rotor, S_rotor{i} ] = jcalc( model.jtype_rotor{i}, q{i}*model.gr{i} );
      S_rotor{i} = S_rotor{i} * model.gr{i};
      vJ_rotor = S_rotor{i} * qd{i};
      Xup_rotor{i} = XJ_rotor * model.Xrotor{i};
      if model.parent(i) == 0
          v_rotor{i} = vJ_rotor;
          c_rotor{i} = zeros(size(a_grav));
      else
          v_rotor{i} = Xup_rotor{i}*v{model.parent(i)} + vJ_rotor;
          c_rotor{i} = crm(v_rotor{i}) * vJ_rotor;
      end
      pA_rotor{i} = crf(v_rotor{i}) * model.I_rotor{i} * v_rotor{i};
  end
end

if nargin == 5
  pA = apply_external_forces( model.parent, Xup, pA, f_ext );
end

for i = model.NB:-1:1
  p = model.parent(i);
  
  % Modified Gauss principle invocation for rotors
  if model.has_rotor(i) 
      U_link{i}  = IA{i} * S{i};
      U_rotor{i} = model.I_rotor{i} * S_rotor{i};
      d{i}  = S{i}' * U_link{i} + S_rotor{i}'*U_rotor{i};

      u{i} = tau{i} - S{i}'*pA{i} - S_rotor{i}'*pA_rotor{i} - U_link{i}'*c{i} - U_rotor{i}'*c_rotor{i};
      U{i} = Xup{i}' * U_link{i} + Xup_rotor{i}' * U_rotor{i};

      if model.parent(i) ~= 0
        Ia = Xup{i}' * IA{i} * Xup{i} + Xup_rotor{i}'*model.I_rotor{i}*Xup_rotor{i};  
        Ia = Ia - U{i}/d{i}*U{i}';

        pa = Xup{i}'*(pA{i} + IA{i}*c{i}) + Xup_rotor{i}'*(pA_rotor{i} + model.I_rotor{i} * c_rotor{i});    
        pa = pa + U{i} * (d{i}\u{i});

        IA{p} = IA{p} + Ia;
        pA{p} = pA{p} + pa;
      end
  else % Usual ABA - slightly modified
      U{i} = IA{i} * S{i};
      d{i} = S{i}' * U{i};
      u{i} = tau{i} - S{i}'*pA{i} - U{i}'*c{i};
      U{i} = Xup{i}'*U{i};
      
      if model.parent(i) ~= 0
        Ia = Xup{i}' * IA{i} * Xup{i};  
        Ia = Ia - U{i}/d{i}*U{i}';

        pa = Xup{i}'*(pA{i} + IA{i}*c{i});    
        pa = pa + U{i} * (d{i}\u{i});

        IA{p} = IA{p} + Ia;
        pA{p} = pA{p} + pa;
      end
  end
end

qdd = zeros(model.NV,1);
for i = 1:model.NB
  if model.parent(i) == 0
     ap= -a_grav; 
  else
     ap = a{model.parent(i)};
  end
  ii = model.vinds{i};
  qdd(ii) = d{i}\(u{i} - U{i}'*ap);
  a{i} = Xup{i} * ap + S{i}*qdd(ii) + c{i};
end