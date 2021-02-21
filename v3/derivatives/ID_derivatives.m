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
if any( model.nv > 1)
    error('ID_derivatives only supports single-DoF joints');
end


a_grav = get_gravity(model);
IC = model.I;
I = model.I;

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i}  = vJ;
    a{i}  = Xup{i}*(-a_grav) + crm(v{i})*vJ + S{i}*qdd{i};
  else
    v{i}  = Xup{i}*v{model.parent(i)} + vJ;
    a{i}  = Xup{i}*a{model.parent(i)} + crm(v{i})*vJ + S{i}*qdd{i};
  end
  Sd{i} = crm(v{i})*S{i};
  Sdd{i}= crm(a{i})*S{i} + crm(v{i})*Sd{i};
  
  BC{i} = factorFunctions(I{i},v{i});
  f{i}  =  I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
end

for i = model.NB:-1:1
  tmp1 = IC{i}*S{i};
  tmp2 = 2*(BC{i}*S{i} + IC{i}*Sd{i});
  tmp3 = 2*BC{i}*Sd{i} + IC{i}*Sdd{i} + crf(S{i})*f{i}; 
  tmp4 = 2*BC{i}'*S{i};
  
  j = i;
  while j > 0
    dtau_dq(j,i)  = S{j}'*tmp3;
    dtau_dq(i,j)  = Sdd{j}'*tmp1+Sd{j}'*tmp4;
    dtau_dqd(j,i) = S{j}'*tmp2;
    dtau_dqd(i,j) = 2*Sd{j}'*tmp1+S{j}'*tmp4;
    
    tmp1 = Xup{j}'*tmp1;
    tmp2 = Xup{j}'*tmp2;
    tmp3 = Xup{j}'*tmp3;
    tmp4 = Xup{j}'*tmp4;
    j  = model.parent(j);
  end
  
  if model.parent(i) > 0
     p = model.parent(i);
     IC{p} = IC{p} + Xup{i}'*IC{i}*Xup{i};
     BC{p} = BC{p} + Xup{i}'*BC{i}*Xup{i};
     f{p}  = f{p}  + Xup{i}'*f{i};
  end
end

