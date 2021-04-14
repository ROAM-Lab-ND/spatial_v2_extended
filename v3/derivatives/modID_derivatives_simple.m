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
  wJ{i} = S{i}*lambda{i};
  
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    ap{i} = zeros(6,1);
    wp{i} = zeros(6,1);
  else
    wp{i} = Xup{i}*w{model.parent(i)};
    ap{i} = Xup{i}*a{model.parent(i)};
  end
  a{i} = ap{i} + S{i}*qdd{i};
  w{i} = wp{i} + wJ{i};
  
  h{i} = I{i}*w{i};
  f{i} = I{i}*a{i};
end

grad_q = repmat( 0*q{1}(1), [size(q,1),size(lambda,2)] );

for i = model.NB:-1:1
   ii = model.vinds{i};
   grad_q(ii,:)  = -S{i}.'*(crf(ap{i})*h{i} + icrf(f{i})*wp{i});
            
   p = model.parent(i);
   if p > 0
      h{p} = h{p}+Xup{i}.'*h{i};
      f{p} = f{p}+Xup{i}.'*f{i};
   end 
end
