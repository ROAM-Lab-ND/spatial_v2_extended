function  [H_derivatives] = H_derivatives( model, q)

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q)
    [q] = confVecToCell(model,q);
end

if sum(model.has_rotor) > 1
    error('H_diff does not support rotors');
end

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  Xup{i} = XJ * model.Xtree{i};
end
IC = model.I;				% composite inertia calculation

H_derivatives = repmat(0*q{1}(1),model.NV,model.NV,model.NV);

for k = model.NB:-1:1
    for k_ind = 1:length(model.vinds{k})
        kk = model.vinds{k}(k_ind);
        sk = S{k}(:,k_ind);
        
        Q = crf(sk)*IC{k} - IC{k}*crm(sk); % Rate of change in IC{k} due to motion of joint k
        Fk = IC{k}*sk; % Other term that shows up in CRBA
        j = k;
        while j > 0
           jj = model.vinds{j};
           F1 = icrf(Fk)*S{j}; % Rate of change in Fk due to motion of joint j
           F2 = Q*S{j};
           i = j;
           while i > 0
               ii = model.vinds{i};
               if i < j
                   H_derivatives(ii,kk,jj) = S{i}.'*F1;    
                   H_derivatives(kk,ii,jj) = S{i}.'*F1;    
               end
              
               if j < k
                H_derivatives(ii,jj,kk) = S{i}.'*F2;    
                H_derivatives(jj,ii,kk) = (S{i}.'*F2).';
               end               
               
               F1 = Xup{i}.'*F1;     F2 = Xup{i}.'*F2;
               i = model.parent(i);
           end
           Fk = Xup{j}'*Fk;     Q = Xup{j}'*Q*Xup{j};
           j = model.parent(j);
        end
    end
    if model.parent(k) > 0
        IC{model.parent(k)} = IC{model.parent(k)} + Xup{k}.'*IC{k}*Xup{k};
    end 
end