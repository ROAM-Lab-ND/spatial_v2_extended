
function Gamma = Christoffel(model,q)

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q) || ~iscell(qd)
    [q] = confVecToCell(model,q);
end
if sum(model.has_rotor) > 1
    error('Christoffel does not support rotors');
end
if any( model.nv > 1)
    error('Christoffel only supports single-DoF joints');
end

IC = model.I;
for i = 1:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
    Xup{i} = XJ * model.Xtree{i};
end

for k = model.NB:-1:1
    B = factorFunctions(IC{k},S{k});
    F = IC{k}*S{k};
    j = k;
    while j > 0
       f1 = B * S{j};
       f2 = B'* S{j};
       f3 = icrf(F)*S{j} - f1;
       i=j;
       while i > 0
           Gamma(i,j,k) = S{i}'*f1;
           Gamma(i,k,j) = Gamma(i,j,k);
           Gamma(j,i,k) = S{i}'*f2;
           Gamma(j,k,i) = Gamma(j,i,k);
           Gamma(k,i,j) = S{i}'*f3;
           Gamma(k,j,i) = Gamma(k,i,j);
           f1 = Xup{i}'*f1;
           f2 = Xup{i}'*f2;
           f3 = Xup{i}'*f3;
           i = model.parent(i);
       end
       B = Xup{j}'*B*Xup{j};
       F = Xup{j}'*F;
       j = model.parent(j);
    end
    if model.parent(k) > 0
        p = model.parent(k);
        IC{p}= IC{p} + Xup{k}'*IC{k}*Xup{k};    
    end
end
