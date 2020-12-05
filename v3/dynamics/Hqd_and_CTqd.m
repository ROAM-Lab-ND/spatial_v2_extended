function  [Hqd, CTqd]  = Hqd_and_CTqd( model, q , qd)

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q)
    [q, qd] = confVecToCell(model,q,qd);
end


for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
  end
  h{i} = model.I{i}*v{i};
  Sd{i} = crm(v{i})*S{i};
  
 
  % Extra data for rotors
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
      Sd_rotor{i} = crm(v_rotor{i})*S_rotor{i};
  end

end

Hqd = zeros(model.NV,1);
CTqd= zeros(model.NV,1);

for i = model.NB:-1:1   
  ii = model.vinds{i};
  Hqd(ii) = S{i}' * h{i};
  CTqd(ii)= Sd{i}'* h{i};
  if model.parent(i) ~= 0
     p = model.parent(i);
     h{p} = h{p} + Xup{i}'*h{i};
  end
  
  if model.has_rotor(i)
      Hqd(ii) = Hqd(ii) + S_rotor{i}' * h_rotor{i};
      CTqd(ii)= CTqd(ii)+ Sd_rotor{i}'* h_rotor{i};
      if model.parent(i) ~= 0
         p = model.parent(i);
         h{p} = h{p} + Xup_rotor{i}'*h_rotor{i};
      end
  end
end

