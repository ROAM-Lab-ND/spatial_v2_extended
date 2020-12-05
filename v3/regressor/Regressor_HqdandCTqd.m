function  [Y_Hqd, Y_CTqd, Y_Hqd_rot, Y_CTqd_rot]  = Regressor_HqdandCTqd( model, q , qd)
% Indirect regressors term  YHqd a = H qd and YCTqd = C^T qd

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q)
    [q, qd] = confVecToCell(model,q,qd);
end

Y_Hqd = zeros(model.NV, 10*model.NB);
Y_Hqd_rot = zeros(model.NV,10*model.NB_rot);

Y_CTqd = zeros(model.NV, 10*model.NB);
Y_CTqd_rot = zeros(model.NV,10*model.NB_rot);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
  end
  Sd{i} = crm(v{i})*S{i};
  hi = individualRegressor(v{i}, 0*v{i});
  param_inds = 10*(i-1) + (1:10);
  j = i;
  while j > 0
    jj = model.vinds{j};
    Y_Hqd(jj, param_inds) = S{j}' * hi;
    Y_CTqd(jj, param_inds)= Sd{j}'* hi;
    hi = Xup{j}' * hi;
	j = model.parent(j);
  end
  
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
    Sd_rotor{i} = crm(v_rotor{i})*S_rotor{i};
    
    hi_rotor = individualRegressor(v_rotor{i}, v_rotor{i}*0);
    
    param_inds = model.rotor_param_inds{i};
    ii = model.vinds{i};
    
    Y_Hqd_rot(ii, param_inds) = S_rotor{i}' * hi_rotor;
    Y_CTqd_rot(ii,param_inds) = Sd_rotor{i}'* hi_rotor;
    
    hi_rotor = Xup_rotor{i}'*hi_rotor;
    j = model.parent(i);
    
    while j > 0
        jj = model.vinds{j};
        Y_Hqd_rot(jj, param_inds) = S{j}' * hi_rotor;
        Y_CTqd_rot(jj,param_inds) = Sd{j}'* hi_rotor;
        hi_rotor = Xup{j}' * hi_rotor;
        j = model.parent(j);
    end
  end
end

