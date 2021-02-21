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

a_grav = get_gravity(model);

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q) || ~iscell(qd)
    [q, qd] = confVecToCell(model,q,qd);
end

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    avp{i} = Xup{i} * -a_grav;
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    avp{i} = Xup{i}*avp{model.parent(i)} + crm(v{i})*vJ;
  end
  fvp{i} = model.I{i}*avp{i} + crf(v{i})*model.I{i}*v{i};
  
  % Extra data for rotors
  if model.has_rotor(i)
      [ XJ_rotor, S_rotor{i} ] = jcalc( model.jtype_rotor{i}, q{i}*model.gr{i} );
      S_rotor{i} = S_rotor{i} * model.gr{i};
      vJ_rotor = S_rotor{i} * qd{i};
      Xup_rotor{i} = XJ_rotor * model.Xrotor{i};
      if model.parent(i) == 0
          v_rotor{i} = vJ_rotor;
          avp_rotor{i} = Xup_rotor{i}*(-a_grav);
      else
          v_rotor{i} = Xup_rotor{i}*v{model.parent(i)} + vJ_rotor;
          avp_rotor{i} = Xup_rotor{i}*avp{model.parent(i)} + crm(v_rotor{i})*vJ_rotor;
      end
      fvp_rotor{i} = model.I_rotor{i}*avp_rotor{i} + crf(v_rotor{i})*model.I_rotor{i}*v_rotor{i};
  end
 
end

if nargin == 4
  fvp = apply_external_forces( model.parent, Xup, fvp, f_ext );
end

C = zeros(model.NV,1);
H = zeros(model.NV);
IC = model.I;				% composite inertia calculation

for i = model.NB:-1:1
  fh = IC{i} * S{i};
  ii = model.vinds{i};
  H(ii,ii) = S{i}' * fh;
  C(ii,1) = S{i}' * fvp{i};
  
  fh = Xup{i}.'* fh;
  
  if model.has_rotor(i)
      H(ii,ii) = H(ii,ii) + S_rotor{i}'* model.I_rotor{i}*S_rotor{i};
      C(ii,1)  = C(ii,1)  + S_rotor{i}'*fvp_rotor{i};
      fh = fh + Xup_rotor{i}'*model.I_rotor{i}*S_rotor{i};
  end
  
  j = i;
  while model.parent(j) > 0
    j = model.parent(j);
    jj = model.vinds{j};
    H(jj,ii) = S{j}' * fh;
    H(ii,jj) = H(jj,ii).';
    fh = Xup{j}.' * fh;
  end
  
  if model.parent(i) ~= 0
    IC{model.parent(i)} = IC{model.parent(i)} + Xup{i}.'*IC{i}*Xup{i};
    fvp{model.parent(i)} = fvp{model.parent(i)} + Xup{i}.'*fvp{i};
    
    if model.has_rotor(i)
        IC{model.parent(i)} = IC{model.parent(i)} + Xup_rotor{i}.'*model.I_rotor{i}*Xup_rotor{i};
        fvp{model.parent(i)} = fvp{model.parent(i)} + Xup_rotor{i}.'* fvp_rotor{i};
    end
  end
end

% info.IC = IC;
% info.fvp = fvp;
% info.avp = avp;
% info.Xup = Xup;
% 
% 
% info.v = v;
% info.Xup_rot = Xup_rot;
% info.v_rot = v_rot;

