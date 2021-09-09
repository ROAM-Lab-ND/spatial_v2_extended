function  [dtau_dq, dtau_dqd] = ID_derivatives( model, q, qd, qdd )

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q) || ~iscell(qd) || ~iscell(qdd)
    [q, qd, qdd] = confVecToCell(model,q,qd,qdd);
end
if sum(model.has_rotor) > 1
    error('ID_derivatives does not support rotors');
end

a_grav = get_gravity(model);
IC = model.I;
I = model.I;

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  Xup{i} = XJ * model.Xtree{i};
  
  if model.parent(i) == 0
    v{i}  = zeros(6,1);
    a{i}  = -a_grav;
    Xup0{i} = Xup{i};
  else
    Xup0{i} = Xup{i}*Xup0{model.parent(i)};
    v{i}  = v{model.parent(i)};
    a{i}  = a{model.parent(i)};
  end
  
  Xdown0{i} = inv(Xup0{i});
  
  S{i} = Xdown0{i}*S{i};
  vJ{i}= S{i}*qd{i};
  aJ{i} = crm(v{i})*vJ{i} + S{i}*qdd{i};
  
  Sd{i} = crm(v{i})*S{i};
  Sdd{i}= crm(a{i})*S{i} + crm(v{i})*Sd{i};
  Sj{i} = 2*Sd{i} + crm(vJ{i})*S{i};
  
  v{i} = v{i} + vJ{i};
  a{i} = a{i} + aJ{i};
  IC{i} = Xup0{i}'*I{i}*Xup0{i};
  
  BC{i} = 2*factorFunctions(IC{i},v{i});
  f{i}  =  IC{i}*a{i} + crf(v{i})*IC{i}*v{i};
end

dtau_dq  = q{1}(1)*0 + zeros(model.NV,model.NV);
dtau_dqd = q{1}(1)*0 + zeros(model.NV,model.NV);

J  = cell2mat(S);
Jd = cell2mat(Sd);
Jdd= cell2mat(Sdd);
Jj = cell2mat(Sj);

for i = model.NB:-1:1
  ii = model.vinds{i};
  
  tmp1(:,ii) = IC{i}*S{i};
  tmp2(:,ii) = BC{i}*S{i} + IC{i}*Sj{i};
  tmp3(:,ii) = BC{i}*Sd{i} + IC{i}*Sdd{i} + icrf(f{i})*S{i}; 
  tmp4(:,ii) = BC{i}.'*S{i};
  
  jj = model.subtree_vinds{i};
  
  dtau_dq(ii,jj)  = J(:,ii).'*tmp3(:,jj);
  dtau_dq(jj,ii)  = tmp1(:,jj).'*Jdd(:,ii)+tmp4(:,jj).'*Jd(:,ii);

  dtau_dqd(ii,jj) = J(:,ii).'*tmp2(:,jj);
  dtau_dqd(jj,ii) =  tmp1(:,jj).'*Jj(:,ii)+tmp4(:,jj).'*J(:,ii);
   
  if model.parent(i) > 0
     p = model.parent(i);
     IC{p} = IC{p} + IC{i};
     BC{p} = BC{p} + BC{i};
     f{p}  = f{p}  + f{i};
  end
end

