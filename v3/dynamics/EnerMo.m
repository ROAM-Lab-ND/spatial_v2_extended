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
    model = model.postProcessModel();
end

[q, qd] = confVecToCell(q,qd);

v = {};

for i = 1:model.NB
  vp = model.getParentVariable(i, v);
  [Xup{i}, S{i}, Sd{i}, v{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i}, qd{i}, vp); 
  
  Ic{i} = model.I{i};
  hc{i} = Ic{i} * v{i};
  KE(i) = 0.5 * v{i}.' * hc{i};
end

ret.Itot = q{1}(1)*0 + zeros(size(Ic{1}));
ret.htot = q{1}(1)*0 + zeros(size(hc{1}));

for i = model.NB:-1:1
  p = model.parent(i);
  
  if p ~= 0
    Ic{p} = Ic{p} + Xup{i}.'*Ic{i}*Xup{i};
    hc{p} = hc{p} + Xup{i}.'*hc{i};
  else
    ret.Itot = ret.Itot + Xup{i}.'*Ic{i}*Xup{i};
    ret.htot = ret.htot + Xup{i}.'*hc{i};
  end
end

a_grav = model.getGravity(model);

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
