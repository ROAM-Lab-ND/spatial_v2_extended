function  [H_diff] = H_diff( model, q)

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q)
    [q] = confVecToCell(model,q);
end

if sum(model.has_rotor) > 1
    error('H_diff does not support rotors');
end

if any( model.nv > 1)
    error('H_diff only supports single-DoF joints');
end



for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  Xup{i} = XJ * model.Xtree{i};
end
IC = model.I;				% composite inertia calculation

H_diff = repmat(0*q{1},model.NV,model.NV,model.NV);

for k = model.NB:-1:1
    Q = crf(S{k})*IC{k} - IC{k}*crm(S{k}); % Rate of change in IC{k} due to motion of joint k
    Fk = IC{k}*S{k}; % Other term that shows up in CRBA
    j = k;
    while j > 0
       F1 = crf(S{j})*Fk; % Rate of change in Fk due to motion of joint j
       F2 = Q*S{j};
       i = j;
       while i > 0
           H_diff(i,k,j) = S{i}'*F1;    H_diff(k,i,j) = S{i}'*F1;    
           H_diff(i,j,k) = S{i}'*F2;    H_diff(j,i,k) = S{i}'*F2;
           
           F1 = Xup{i}'*F1;     F2 = Xup{i}'*F2;
           i = model.parent(i);
       end
       Fk = Xup{j}'*Fk;     Q = Xup{j}'*Q*Xup{j};
       j = model.parent(j);
    end
    if model.parent(k) > 0
        IC{model.parent(k)} = IC{model.parent(k)} + Xup{k}'*IC{k}*Xup{k};
    end
end