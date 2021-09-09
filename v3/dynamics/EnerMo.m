function  ret = EnerMo( model, q, qd )

% EnerMo  calculate energy, momentum and related quantities
% EnerMo(robot,q,qd)  returns a structure containing the fields KE, PE,
% htot, Itot, mass, cm and vcm.  These fields contain the kinetic and
% potential energies of the whole system, the total spatial momentum, the
% total spatial inertia, total mass, position of centre of mass, and the
% linear velocity of centre of mass, respectively.  Vector quantities are
% expressed in base coordinates.  PE is defined to be zero when cm is
% zero.

if ~isfield(model,'nq')
    model = postProcessModel(model);
end

[q, qd] = confVecToCell(model,q,qd);

% if sum(model.has_rotor) > 1
%     error('EnerMo does not support rotors');
% end


for i = 1:model.NB
  [ XJ, S ] = jcalc( model.jtype{i}, q{i} );
  vJ = S*qd{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
  end
  
  Ic{i} = model.I{i};
  hc{i} = Ic{i} * v{i};
  KE(i) = 0.5 * v{i}.' * hc{i};
  if model.has_rotor(i)
      [ XJ_rotor, S_rotor{i} ] = jcalc( model.jtype_rotor{i}, q{i}*model.gr{i} );
      S_rotor{i} = S_rotor{i} * model.gr{i};
      vJ_rotor = S_rotor{i} * qd{i};
      Xup_rotor{i} = XJ_rotor * model.Xrotor{i};
      if model.parent(i) == 0
          v_rotor{i} = vJ_rotor;
      else
          v_rotor{i} = Xup_rotor{i}*v{model.parent(i)} + vJ_rotor;
      end
      h_rotor{i} = model.I_rotor{i}*v_rotor{i};
      KE(i) = KE(i) + 1/2*h_rotor{i}.'*v_rotor{i};
  end
end

ret.Itot = q{1}(1)*0 + zeros(size(Ic{1}));
ret.htot = q{1}(1)*0 + zeros(size(hc{1}));

for i = model.NB:-1:1
  if model.parent(i) ~= 0
    Ic{model.parent(i)} = Ic{model.parent(i)} + Xup{i}.'*Ic{i}*Xup{i};
    hc{model.parent(i)} = hc{model.parent(i)} + Xup{i}.'*hc{i};
    
      if model.has_rotor(i) % Modified backward pass
          Ic{model.parent(i)} = Ic{model.parent(i)} + Xup_rotor{i}.'*model.I_rotor{i}*Xup_rotor{i};
          hc{model.parent(i)} = hc{model.parent(i)} + Xup_rotor{i}.'* h_rotor{i};
      end
  else
    ret.Itot = ret.Itot + Xup{i}.'*Ic{i}*Xup{i};
    ret.htot = ret.htot + Xup{i}.'*hc{i};
    if model.has_rotor(i) % Modified backward pass
        ret.Itot = ret.Itot + Xup_rotor{i}.'*model.I_rotor{i}*Xup_rotor{i};
        ret.htot = ret.htot + Xup_rotor{i}.'* h_rotor{i};
    end
  end
end

a_grav = get_gravity(model);

if length(a_grav) == 6
  g = a_grav(4:6);			% 3D linear gravitational accn
  h = ret.htot(4:6);			% 3D linear momentum
else
  g = a_grav(2:3);			% 2D gravity
  h = ret.htot(2:3);			% 2D linear momentum
end

[mass, cm] = mcI(ret.Itot);

ret.KE = sum(KE);
ret.PE = - mass * dot(cm,g);
ret.mass = mass;
ret.cm = cm;
ret.vcm = h / mass;
