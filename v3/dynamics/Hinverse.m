function  [Hinv] = Hinverse( model, q)
%
% Factored ABA algorithm for computing the inverse of the mass matrix H
%
% Output Hinv satisfies: Hinv = inv(H)
%
% This condition is equivalent to:
%               q_ddot = Hinv * tau when q_dot=0 and a_grav = 0
% 

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q)
    [q] = confVecToCell(model,q);
end
if sum(model.has_rotor) > 1
    error('Hinverse does not support rotors');
end


IA = model.I; % Cell array for Artiulated-Body Inertias
F = repmat({zeros(6,model.NV)},model.NV,1); % Satisfies F{i}*tau = pA{i} from usual ABA
P = repmat({zeros(6,model.NV)},model.NV,1); % Satisfies P{i}*tau = a{i}  from usual ABA

Hinv = zeros(model.NV); % Inerse of Mass Matrix

% Outward Pass
for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  Xup{i} = XJ * model.Xtree{i};
end

% Inward Pass
for i = model.NB:-1:1   
   ii = model.vinds{i};
   U{i} = IA{i}*S{i};
   D{i} = S{i}'*U{i};
   Hinv(ii,:) = -D{i}\S{i}'*F{i};  
   Hinv(ii,ii) = Hinv(i,i) + inv(D{i});
   % Note: At this point, Hinv(i,:) is not equal to its final value.
   %       However, it does satisfy Hinv(i,:)*tau = D{i}\u{i} from usual ABA
   
   p = model.parent(i);
   if p > 0   
      Fa = F{i} + U{i}*Hinv(ii,:);     % Bias force transmitted to predecessor
      Ia = IA{i} - U{i}*(D{i}\U{i}'); % Articulated inertia transmitted to predeceesor
      
      F{p} = F{p} + Xup{i}'*Fa;
      IA{p} = IA{p} + Xup{i}'*Ia*Xup{i};
   end
end

% Outward Pass
for i = 1:model.NB
   p = model.parent(i);
   ii = model.vinds{i};
   if p > 0
       Hinv(ii,:) = Hinv(ii,:) - D{i}\U{i}'*Xup{i}*P{p};
       P{i} = Xup{i}*P{p} + S{i}*Hinv(ii,:); 
   else
       P{i} = S{i}*Hinv(ii,:);
   end
end