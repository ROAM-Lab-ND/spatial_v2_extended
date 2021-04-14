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
  v{i} = vp{i} + vJ{i};
  a{i} = ap{i} + crm(v{i})*vJ{i} + S{i}*qdd{i};
  w{i} = wp{i} + wJ{i};
  
  %out = out + w{i}.'*(I{i}*a{i} + crf(v{i})*I{i}*v{i});
  
  h{i} = I{i}*w{i};
  %z{i} = (-I{i}*crm(v{i}) - icrf(I{i}*v{i}))*w{i};
  z{i} = -( crf(w{i})*I{i} - I{i}*crm(w{i}))*v{i};
  f{i} = I{i}*a{i} + crf(v{i})*I{i}*v{i};
  
end

grad_q = repmat( 0*q{1}(1), [size(qd,1),size(lambda,2)] );
grad_qd = repmat(0*0*q{1}(1),[size(qd,1),size(lambda,2)] );

for i = model.NB:-1:1
   ii = model.vinds{i};
   z{i} = z{i} + crf(vJ{i})*h{i};    
   grad_qd(ii,:) = S{i}.'*(z{i}-crf(v{i})*h{i});
   grad_q(ii,:)  = -S{i}.'*(crf(vp{i})*z{i} +crf(ap{i})*h{i} + icrf(f{i})*wp{i});
   p = model.parent(i);
   if p > 0
      z{p} = z{p}+Xup{i}.'*z{i};
      h{p} = h{p}+Xup{i}.'*h{i};
      f{p} = f{p}+Xup{i}.'*f{i};
   end 
end
total_grad = [ grad_q ; grad_qd];
