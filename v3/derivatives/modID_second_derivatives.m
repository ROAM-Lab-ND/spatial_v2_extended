function  [grad_q, grad_qd, H_qq, H_qdqd, H_qdq, dtau_dq, dtau_dqd] = modID_second_derivatives( model, q, qd, qdd, lambda)

%%% Forward mode chain rule applied to the reverse mode accmulation in 
% mod_ID_derivatives
% 
% It's not pretty. Must be something more elegant to find...
%

if ~isfield(model,'nq')
    model = postProcessModel(model);
end

if sum(model.has_rotor) > 1
    error('modID does not support rotors');
end

a_grav = get_gravity(model);
if ~iscell(q)
    [q, qd, qdd, lambda] = confVecToCell(model,q,qd,qdd, lambda);
end

%out = zeros(size(lambda,2),1);

dv_dq_p = repmat({zeros(6,model.NB)},model.NB,1);
dv_dqd_p = repmat({zeros(6,model.NB)},model.NB,1);

da_dq_p = repmat({zeros(6,model.NB)},model.NB,1);
da_dqd_p = repmat({zeros(6,model.NB)},model.NB,1);

wp=repmat({zeros(6,1)},model.NB,1);
dw_dq_p = repmat({zeros(6,model.NB)},model.NB,1);


I = model.I;

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ{i} = S{i}*qd{i};
  wJ = S{i}*lambda{i};
  
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    ap{i} = -a_grav;
    vp{i} = zeros(6,1);
    wp{i} = zeros(6,1);
  else
    vp{i} = Xup{i}*v{model.parent(i)};
    wp{i} = Xup{i}*w{model.parent(i)};
    ap{i} = Xup{i}*a{model.parent(i)};
    
    dv_dq_p{i} = Xup{i}*dv_dq{model.parent(i)} ;
    dv_dq_p{i}(:,i) = -crm(S{i})*vp{i};
    
    
    da_dq_p{i} = Xup{i}*da_dq{model.parent(i)} ;
    da_dq_p{i}(:,i) = -crm(S{i})*ap{i};
    
    dw_dq_p{i} = Xup{i}*dw_dq{model.parent(i)} ;
    dw_dq_p{i}(:,i) = -crm(S{i})*wp{i};
    
    
    dv_dqd_p{i} =  Xup{i}*dv_dqd{model.parent(i)} ;
    da_dqd_p{i} =  Xup{i}*da_dqd{model.parent(i)} ;
    
  end
  v{i} = vp{i} + vJ{i};
  a{i} = ap{i} + crm(v{i})*vJ{i} + S{i}*qdd{i};
  w{i} = wp{i} + wJ;
  
  dv_dq{i} = dv_dq_p{i};
  
  da_dq{i} = da_dq_p{i};
  da_dq{i} = da_dq{i} - crm(vJ{i})*dv_dq{i};
  
  %
  dw_dq{i} = dw_dq_p{i};
  
  dv_dqd{i} = dv_dqd_p{i};
  dv_dqd{i}(:,i) = S{i};
  
  da_dqd{i} = da_dqd_p{i}-crm(vJ{i})*dv_dqd{i};
  da_dqd{i}(:,i) = crm(v{i})*S{i};
  
  % out = out + w{i}'*(I{i}*a{i} + crf(v{i})*I{i}*v{i});
  
  %
  h{i} = I{i}*w{i};
  dh_dq{i} = I{i}*dw_dq{i};
  
  %
  z{i} = I{i}*crm(w{i})*v{i} - crf(w{i})*I{i}*v{i};
  
  dz_dq{i} = I{i}*crm(w{i})*dv_dq{i} - I{i}*crm(v{i})*dw_dq{i} ...
            - crf(w{i})*I{i}*dv_dq{i} - icrf(I{i}*v{i})*dw_dq{i};
        
  dz_dqd{i} = I{i}*crm(w{i})*dv_dqd{i} - crf(w{i})*I{i}*dv_dqd{i};
        
  f{i} = I{i}*a{i} + crf(v{i})*I{i}*v{i};
  df_dq{i} = I{i}*da_dq{i}  + crf(v{i})*I{i}*dv_dq{i} +icrf(I{i}*v{i})*dv_dq{i};
  df_dqd{i}= I{i}*da_dqd{i} + crf(v{i})*I{i}*dv_dqd{i}+icrf(I{i}*v{i})*dv_dqd{i};
  
end

for i = model.NB:-1:1
   z{i} = z{i} + crf(vJ{i})*h{i};
   
   dz_dq{i} = dz_dq{i} + crf(vJ{i})*dh_dq{i};
   dz_dqd{i}(:,i) = dz_dqd{i}(:,i) + crf(S{i})*h{i};
    
   grad_qd(i) = S{i}'*(z{i}-crf(v{i})*h{i});
   grad_q(i)  = -S{i}'*(crf(vp{i})*z{i} +crf(ap{i})*h{i} + crf(wp{i})*f{i});
   
   H_qdq(i,:) = S{i}'*(dz_dq{i}-crf(v{i})*dh_dq{i}-icrf(h{i})*dv_dq{i});
   
   H_qdqd(i,:) = S{i}'*(dz_dqd{i} -icrf(h{i})*dv_dqd{i});
   
   H_qq(i,:) = -S{i}'*(crf(vp{i})*dz_dq{i} + crf(ap{i})*dh_dq{i} + crf(wp{i})*df_dq{i} + ...
                        icrf(z{i})*dv_dq_p{i}+ icrf(h{i})*da_dq_p{i} + icrf(f{i})*dw_dq_p{i});                 
   
   dtau_dq(i,:) = S{i}'*df_dq{i};
   dtau_dqd(i,:) = S{i}'*df_dqd{i};
                    
   p = model.parent(i);
   if p > 0
      z{p} = z{p}+Xup{i}'*z{i};
      dz_dq{p} = dz_dq{p} + Xup{i}'*dz_dq{i};
      dz_dq{p}(:,i) = dz_dq{p}(:,i) - (crm(S{i})*Xup{i})'*z{i};
      
      dz_dqd{p}= dz_dqd{p} + Xup{i}'*dz_dqd{i};
      
      h{p} = h{p}+Xup{i}'*h{i};
      dh_dq{p} = dh_dq{p} + Xup{i}'*dh_dq{i};
      dh_dq{p}(:,i) = dh_dq{p}(:,i) - (crm(S{i})*Xup{i})'*h{i};
      
      f{p} = f{p}+Xup{i}'*f{i};
      df_dq{p} = df_dq{p} + Xup{i}'*df_dq{i};
      df_dq{p}(:,i) = df_dq{p}(:,i) - (crm(S{i})*Xup{i})'*f{i};
      
      df_dqd{p}= df_dqd{p}+ Xup{i}'*df_dqd{i};
   end 
end