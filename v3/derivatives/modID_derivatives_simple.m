function  [grad_q] = modID_derivatives_simple( model, q, qdd, lambda )
%%% Simplified version of modID_derivatives when gravity off, and qd=0

if ~isfield(model,'nq')
    model = postProcessModel(model);
end

if sum(model.has_rotor) > 1
    error('modID does not support rotors');
end

if ~iscell(q)
    [q, qdd, lambda] = confVecToCell(model,q,qdd, lambda);
end

% Calculate the derivatives using reverse mode chain rule
for i = 1:model.NB
  I{i} = q{1}(1)*0+ model.I{i};
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    X0{i} = Xup{i};
    a{i} = zeros(6,1);
    w{i} = zeros(6,1);
  else
    X0{i} = Xup{i}*X0{model.parent(i)};
    w{i} = w{model.parent(i)};
    a{i} = a{model.parent(i)};
  end
  wp{i} = w{i};
  
  S{i}   = X0{i}\S{i};
  I{i}   = X0{i}.'*I{i}*X0{i};
  Ydd{i} = crm(a{i})*S{i};
  
  a{i} = a{i} + S{i}*qdd{i};
  w{i} = w{i} + S{i}*lambda{i};
  
  h{i} = I{i}*w{i};
  f{i} = I{i}*a{i};
end

grad_q = repmat( 0*q{1}(1), [size(q,1),size(lambda,2)] );

for i = model.NB:-1:1
   ii = model.vinds{i};
   grad_q(ii,:)  = Ydd{i}.'*h{i} + (icrf(f{i})*S{i})'*wp{i};
            
   p = model.parent(i);
   if p > 0
      h{p} = h{p}+h{i};
      f{p} = f{p}+f{i};
   end 
end
