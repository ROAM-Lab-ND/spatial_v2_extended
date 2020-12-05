function  [tau,info] = ID_SlotineLI_rotoror( model, q, qd , qd_r, qdd, factorFunction)

if nargin == 5
    factorFunction = @(I,v)(factorFunctions(I,v));
end

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q)
    [q,qd,qd_r,qdd] = confVecToCell(model,q,qd,qd_r,qdd);
    
end

a_grav = get_gravity(model);
p = model.parent;

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  wJ = S{i}*qd_r{i};
  Xup{i} = XJ * model.Xtree{i};
  if p(i) == 0
    v{i} = vJ;
    w{i} = wJ;
    wd{i}= Xup{i}*(-a_grav) + S{i} * qdd{i} + crm(v{i})*wJ;
  else
    v{i} = Xup{i}*v{p(i)} + vJ;
    w{i} = Xup{i}*w{p(i)} + wJ;
    wd{i}= Xup{i}*wd{p(i)} + S{i} * qdd{i} + crm(v{i})*wJ;
  end
  f{i} = model.I{i}*wd{i} + factorFunction(model.I{i}, v{i})*w{i};
  
  % Extra data for rotors
  if model.has_rotor(i)
      [ XJ_rotor, S_rotor{i} ] = jcalc( model.jtype_rotor{i}, q{i}*model.gr{i} );
      S_rotor{i} = S_rotor{i} * model.gr{i};
      vJ_rotor = S_rotor{i} * qd{i};
      wJ_rotor = S_rotor{i} * qd_r{i};
       
      Xup_rotor{i} = XJ_rotor * model.Xrotor{i};
      if p(i) == 0
          v_rotor{i} = vJ_rotor;
          w_rotor{i} = wJ_rotor;
          wd_rotor{i} = Xup_rotor{i}*(-a_grav) + S_rotor{i}*qdd{i};
      else
          v_rotor{i}  = Xup_rotor{i}*v{p(i)} + vJ_rotor;
          w_rotor{i}  = Xup_rotor{i}*w{p(i)} + wJ_rotor;
          wd_rotor{i} = Xup_rotor{i}*wd{p(i)} + S_rotor{i}*qdd{i} + crm(v_rotor{i})*wJ_rotor;
      end
      f_rotor{i} = model.I_rotor{i}*wd_rotor{i} + factorFunction(model.I_rotor{i}, v_rotor{i})*w_rotor{i};
  end
 
end

tau = zeros(model.NV,1);
for i = model.NB:-1:1
  ii = model.vinds{i};
  tau(ii) = S{i}' * f{i} ;
  
  if p(i) ~= 0
    f{p(i)} = f{p(i)} + Xup{i}'*f{i} ;
  end
      
  if model.has_rotor(i) % Modified backward pass
      tau(ii) = tau(ii) + S_rotor{i}' * f_rotor{i};
      if p(i) ~= 0
        f{p(i)} = f{p(i)} +  Xup_rotor{i}'*f_rotor{i};
      end   
  end
end
