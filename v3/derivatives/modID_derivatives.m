function  [grad_q, grad_qd, total_grad] = modID_derivatives( model, q, qd, qdd, lambda )

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

% Calculate the derivatives using reverse mode chain rule
for i = 1:model.NB
  I{i} = q{1}(1)*0+ model.I{i};
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ{i} = S{i}*qd{i};
  wJ{i} = S{i}*lambda{i};
  
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    vp{i} = zeros(6,1);
    ap{i} = Xup{i}*(-a_grav);
    wp{i} = zeros(6,1);
  else
    vp{i} = Xup{i}*v{model.parent(i)};
    wp{i} = Xup{i}*w{model.parent(i)};
    ap{i} = Xup{i}*a{model.parent(i)};
  end
  
  Yd{i}   = crm(vp{i})*S{i};
  Ydd{i}  = crm(vp{i})*Yd{i} + crm(ap{i})*S{i};
  Psid{i} = 2*Yd{i}+crm(vJ{i})*S{i};

  
  v{i} = vp{i} + vJ{i};
  a{i} = ap{i} + crm(v{i})*vJ{i} + S{i}*qdd{i};
  f{i} = I{i}*a{i} + crf(v{i})*I{i}*v{i};
  
  w{i} = wp{i} + wJ{i};
  h{i} = I{i}*w{i}; 
  B{i} = factorFunctions(I{i},v{i});
  z{i} = B{i}.'*w{i};
end

grad_q = repmat( 0*q{1}(1), [size(qd,1),size(lambda,2)] );
grad_qd = repmat(0*0*q{1}(1),[size(qd,1),size(lambda,2)] );

for i = model.NB:-1:1
   ii = model.vinds{i};    
   grad_qd(ii,:) = Psid{i}.'*h{i} + 2*S{i}.'*z{i};
   grad_q(ii,:)  = (icrf(f{i})*S{i}).'*wp{i} + 2*Yd{i}.'*z{i} + Ydd{i}.'*h{i};
   p = model.parent(i);
   if p > 0
      z{p} = z{p}+Xup{i}.'*z{i};
      h{p} = h{p}+Xup{i}.'*h{i};
      f{p} = f{p}+Xup{i}.'*f{i};
   end 
end
total_grad = [ grad_q ; grad_qd];
