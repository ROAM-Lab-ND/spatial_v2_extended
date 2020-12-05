function  [Y, Y_rot] = RegressorClassical( model, q, qd,qdd)

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q)
    [q,qd,qdd] = confVecToCell(model,q,qd,qdd);
    
end

Y = zeros(model.NV, model.NB*10);
Y_rot = zeros(model.NV, model.NB_rot*10);

p = model.parent;
a_grav = get_gravity(model);
for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  Xup{i} = XJ * model.Xtree{i};
  if p(i) == 0
    v{i} = vJ;
    a{i}= -Xup{i}*a_grav + S{i} * qdd{i} + crm(v{i})*vJ;
  else
    v{i} = Xup{i}*v{p(i)} + vJ;
    a{i}= Xup{i}*a{p(i)} + S{i} * qdd{i} + crm(v{i})*vJ;
  end
  F{i} = individualRegressor(a{i}, v{i});
  
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
      F_rotor{i} = individualRegressor(a_rotor{i}, v_rotor{i});
  end  
end



for i = model.NB:-1:1
  ii = model.vinds{i};
  param_inds = 10*(i-1)+1 : 10*i;
  Y(ii,param_inds) = S{i}' * F{i};
  Fi = Xup{i}' * F{i};
  j = i;
  j = model.parent(j);
  while j ~= 0
      jj = model.vinds{j};
      Y(jj,param_inds) = S{j}' * Fi;
      Fi = Xup{j}' * Fi;
      j = model.parent(j);
  end
  
  if model.has_rotor(i)
      param_inds = model.rotor_param_inds{i};
      
      Y_rot(ii,param_inds) = S_rotor{i}' * F_rotor{i};
      Fi = Xup_rotor{i}' * F_rotor{i};
      j = i;
      j = model.parent(j);
      while j ~= 0
          jj = model.vinds{j};
          Y_rot(jj,param_inds) = S{j}' * Fi;
          Fi = Xup{j}' * Fi;
          j = model.parent(j);
      end
  end
end

