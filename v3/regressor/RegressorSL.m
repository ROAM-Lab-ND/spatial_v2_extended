function  [Y, Y_rot] = RegressorSL_rotor( model, q, qd,qd_r,qdd, factorFunction)

if nargin == 5
    factorFunction = @(I,v)(factorFunctions(I,v));
end
% Xup = repmat({zeros(6,6)},model.NB,1);
% Xrot = repmat({zeros(6,6)},model.NB,1);
% X0 = repmat({zeros(6,6)},model.NB,1);
% X0rot = repmat({zeros(6,6)},model.NB,1);
% 
% v   = repmat({zeros(6,1)},model.NB,1);
% w   = repmat({zeros(6,1)},model.NB,1);
% wd  = repmat({zeros(6,1)},model.NB,1);
% vrot   = repmat({zeros(6,1)},model.NB,1);
% wrot   = repmat({zeros(6,1)},model.NB,1);
% wdrot  = repmat({zeros(6,1)},model.NB,1);
% F  = repmat({zeros(6,10)},model.NB,1);
% Frot  = repmat({zeros(6,10)},model.NB,1);

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q)
    [q,qd,qd_r,qdd] = confVecToCell(model,q,qd,qd_r,qdd);
    
end

Y = zeros(model.NV, model.NB*10);
Y_rot = zeros(model.NV, model.NB_rot*10);

p = model.parent;
a_grav = get_gravity(model);
for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  wJ= S{i}*qd_r{i};
  Xup{i} = XJ * model.Xtree{i};
  
  if p(i) == 0
    v{i} = vJ;
    w{i} = wJ;
    wd{i}= -Xup{i}*a_grav + S{i} * qdd{i} + crm(v{i})*wJ;
  else
    v{i} = Xup{i}*v{p(i)} + vJ;
    w{i} = Xup{i}*w{p(i)} + wJ;
    wd{i}= Xup{i}*wd{p(i)} + S{i} * qdd{i} + crm(v{i})*wJ;
  end
  for k = 1:10
     ak = zeros(10,1); ak(k) = 1;
     Ik = inertiaVecToMat(ak);
     F{i}(:,k)    = Ik * wd{i}    + factorFunction(Ik, v{i})    * w{i};
  end
  
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
      
      for k = 1:10
         ak = zeros(10,1); ak(k) = 1;
         Ik = inertiaVecToMat(ak);
         F_rotor{i}(:,k) = Ik * wd_rotor{i} + factorFunction(Ik, v_rotor{i}) * w_rotor{i};
      end
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

