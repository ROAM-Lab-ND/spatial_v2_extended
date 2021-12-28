function  [q_dot, p_dot]  = HamiltonianDynamics( model, q , p, tau, f_ext)

if ~isfield(model,'nq')
    model = postProcessModel(model);
end

model_no_grav = model;
model_no_grav.gravity = [0 0 0]';
qd    = FDab( model_no_grav , q, zeros(model.NV,1), p);
if nargin == 5
    tau_g = ID(model, q, zeros(model.NV,1), zeros(model.NV,1), f_ext);
else
    tau_g = ID(model, q, zeros(model.NV,1), zeros(model.NV,1));
end

if ~iscell(q)
    [q, qd, p] = confVecToCell(model,q,qd,p);
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

CTqd= q{1}(1)*0 + zeros(model.NV,1);
q_dot = q{1}(1)*0 + zeros(model.NV,1);

for i = model.NB:-1:1   
  ii = model.vinds{i};
  q_dot(ii) = qd{i};
  CTqd(ii)= Sd{i}'* h{i};
  if model.parent(i) ~= 0
     p = model.parent(i);
     h{p} = h{p} + Xup{i}'*h{i};
  end
  if model.has_rotor(i)
      CTqd(ii)= CTqd(ii)+ Sd_rotor{i}'* h_rotor{i};
      if model.parent(i) ~= 0
         p = model.parent(i);
         h{p} = h{p} + Xup_rotor{i}'*h_rotor{i};
      end
  end
end

p_dot = CTqd + tau - tau_g;
