function  derivs = modID_second_derivatives( model, q, qd, qdd, lambda)

%%% Forward mode chain rule applied to the reverse mode accmulation in 
% mod_ID_derivatives. 
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

dv_dq_p = repmat({zeros(6,model.NV)},model.NB,1);
dv_dqd_p = repmat({zeros(6,model.NV)},model.NB,1);

da_dq_p = repmat({zeros(6,model.NV)},model.NB,1);
da_dqd_p = repmat({zeros(6,model.NV)},model.NB,1);

wp=repmat({zeros(6,1)},model.NB,1);
dw_dq_p = repmat({zeros(6,model.NV)},model.NB,1);


I = model.I;

for i = 1:model.NB
  
  ii = model.vinds{i};
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ{i} = S{i}*qd{i};
  wJ = S{i}*lambda{i};
  
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    ap{i} = Xup{i}*-a_grav;
    vp{i} = zeros(6,1);
    wp{i} = zeros(6,1);
    
    da_dq_p{i} = zeros(6,model.NV);
    da_dq_p{i}(:,ii) = crm(ap{i})*S{i};
    
  else
    vp{i} = Xup{i}*v{model.parent(i)};
    wp{i} = Xup{i}*w{model.parent(i)};
    ap{i} = Xup{i}*a{model.parent(i)};
    
    dv_dq_p{i} = Xup{i}*dv_dq{model.parent(i)} ;
    dv_dq_p{i}(:,ii) = crm(vp{i})*S{i};
    
    
    da_dq_p{i} = Xup{i}*da_dq{model.parent(i)} ;
    da_dq_p{i}(:,ii) = crm(ap{i})*S{i};
    
    dw_dq_p{i} = Xup{i}*dw_dq{model.parent(i)} ;
    dw_dq_p{i}(:,ii) = crm(wp{i})*S{i};
    
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
  dv_dqd{i}(:,ii) = S{i};
  
  da_dqd{i} = da_dqd_p{i}-crm(vJ{i})*dv_dqd{i};
  da_dqd{i}(:,ii) = da_dqd{i}(:,ii) + crm(v{i})*S{i};
  
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

grad_q = zeros(model.NV,1);
grad_qd = zeros(model.NV,1);

for i = model.NB:-1:1
   ii = model.vinds{i};
   z{i} = z{i} + crf(vJ{i})*h{i};
   
   dz_dq{i} = dz_dq{i} + crf(vJ{i})*dh_dq{i};
   dz_dqd{i}(:,ii) = dz_dqd{i}(:,ii) + icrf(h{i})*S{i};
    
   grad_qd(ii) = S{i}'*(z{i}-crf(v{i})*h{i});
   grad_q(ii)  = -S{i}'*(crf(vp{i})*z{i} +crf(ap{i})*h{i} + crf(wp{i})*f{i});
   
   H_qdq(ii,:) = S{i}'*(dz_dq{i}-crf(v{i})*dh_dq{i}-icrf(h{i})*dv_dq{i});
   
   H_qdqd(ii,:) = S{i}'*(dz_dqd{i} -icrf(h{i})*dv_dqd{i});
   
   H_qq(ii,:) = -S{i}'*(crf(vp{i})*dz_dq{i} + crf(ap{i})*dh_dq{i} + crf(wp{i})*df_dq{i} + ...
                        icrf(z{i})*dv_dq_p{i}+ icrf(h{i})*da_dq_p{i} + icrf(f{i})*dw_dq_p{i});                 
   
   dtau_dq(ii,:) = S{i}'*df_dq{i};
   dtau_dqd(ii,:) = S{i}'*df_dqd{i};
                    
   p = model.parent(i);
   if p > 0
      z{p} = z{p}+Xup{i}'*z{i};
      dz_dq{p} = dz_dq{p} + Xup{i}'*dz_dq{i};
      dz_dq{p}(:,ii) = dz_dq{p}(:,ii) + Xup{i}'*icrf(z{i})*S{i};
      
      dz_dqd{p}= dz_dqd{p} + Xup{i}'*dz_dqd{i};
      
      h{p} = h{p}+Xup{i}'*h{i};
      dh_dq{p} = dh_dq{p} + Xup{i}'*dh_dq{i};
      dh_dq{p}(:,ii) = dh_dq{p}(:,ii) + Xup{i}'*icrf(h{i})*S{i};
      
      f{p} = f{p}+Xup{i}'*f{i};
      df_dq{p} = df_dq{p} + Xup{i}'*df_dq{i};
      df_dq{p}(:,ii) = df_dq{p}(:,ii) + Xup{i}'*icrf(f{i})*S{i};
      
      df_dqd{p}= df_dqd{p}+ Xup{i}'*df_dqd{i};
   end 
end

derivs.dtau_dq  = dtau_dq ;
derivs.dtau_dv  = dtau_dqd;
derivs.dmod_dq  = grad_q ;
derivs.dmod_dv  = grad_qd;
derivs.dmod_dqq = H_qq;
derivs.dmod_dvv = H_qdqd;
derivs.dmod_dqv = H_qdq.';
