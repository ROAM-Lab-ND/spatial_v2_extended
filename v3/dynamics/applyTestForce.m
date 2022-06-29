function [qdd, lambda_inv, qdd_from_subqdd, D] = applyTestForce(model,q, body_number,test_force)
    

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q) || ~iscell(qd)
    [q] = confVecToCell(model,q);
end

% Kinematics
for i = 1:model.NB
  [Xup{i}, S{i}] = model.joint{i}.kinematics(model.Xtree{i},q{i});
  IA{i} = model.I{i};
end

% ABA Step
for i = model.NB:-1:1
  p = model.parent(i);
  ii = model.vinds{i};
  
  U{i}  = IA{i} * S{i};
  d{i}  = S{i}.' * U{i};
  
  D(ii,ii) = d{i};
  
  U{i} = Xup{i}.' * U{i};

  % Force propagators
  Chiup{i} =  (Xup{i} - (S{i}/d{i})*U{i}');
  
  if p ~= 0
    Ia = Xup{i}.' * IA{i} * Xup{i} - U{i}/d{i}*U{i}.';
    IA{p} = IA{p} + Ia;
  end
end

qdd_from_subqdd = zeros(model.NV);
% Update tmp quantities - only needs done once so could be precomputed.
for i = 1:model.NB
  ii = model.vinds{i};  
  qdd_from_subqdd(ii,ii) = eye( length(ii) );
  j = model.parent(i);
  
  % Here we are constructing the basis matrix Psi for the torques
  % Such that f_{i} = Psi{i} * tau{i} + constraints 
  % We set it up so that S{i}'*Psi{i} = Identity. 
  % [Note, this matrix is only unique up to any additive matrix whose range 
  % is orthogonal to the range of S.]
  % Instead of doing a QR, we should probaly bake this into the joint model
  % class down the road so that it is more physically meaningful.
  dofs = length(ii);
  [Q, ~] = qr(S{i}); % Use a QR to get a basis for constrained modes.
  
  S_free_and_constrained = [S{i} Q(:,dofs+1:end)];
  Psi{i} = inv(S_free_and_constrained)';
  Psi{i} = Psi{i}(:,1:dofs);
       
  F = (Chiup{i}'-Xup{i}')*Psi{i};
  while j > 0
    jj = model.vinds{j};
    qdd_from_subqdd(ii,jj) = F.'*S{j};
    F = Chiup{j}'*F;
    j = model.parent(j);
  end
end

F = test_force;
lambda_inv = 0;
i = body_number;
qdd = zeros(model.NV,1);
while i > 0
    ii = model.vinds{i};
    tmp = S{i}'*F;
    lambda_inv = lambda_inv + tmp.'*(d{i}\tmp);
    qdd = qdd + qdd_from_subqdd(:,ii)*(d{i}\tmp);
    F = Chiup{i}'*F;
    i = model.parent(i);
end