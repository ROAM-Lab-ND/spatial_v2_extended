function  [C, Hdot, H ] = CoriolisMatrix( model, q, qd, factorFunction)

if nargin == 3
    factorFunction = @(I,v)(factorFunctions(I,v));
end

if ~isfield(model,'nq')
    model = model.postProcessModel();
end
[q, qd] = model.confVecToCell(q,qd);

I = model.I;

v = {};
for i = 1:model.NB
    vp = model.getParentVariable(i, v);
    [Xup{i}, S{i}, Sdot{i}, v{i}] = model.joint{i}.kinematics(model.Xtree{i}, q{i}, qd{i}, vp);
    BC{i} = factorFunction(I{i}, v{i});    
end

for i =1:model.NB
    IC{i} = q{1}(1)*0 + model.I{i};
end

C    = q{1}(1)*0 + zeros(model.NV);
H    = q{1}(1)*0 + zeros(model.NV);
Hdot = q{1}(1)*0 + zeros(model.NV);

for j = model.NB:-1:1
    jj = model.vinds{j};
   
    f1 = IC{j} * Sdot{j} + BC{j} * S{j};
    f2 = IC{j} * S{j};
    f3 = BC{j}.' * S{j};
    
    C(jj,jj) = S{j}.' * f1;
    H(jj,jj) = S{j}.' * f2;
    Hdot(jj,jj) = Sdot{j}.'*f2 + S{j}.'*(f1 + f3);
    
    f1 = Xup{j}.'*f1;
    f2 = Xup{j}.'*f2;
    f3 = Xup{j}.'*f3;

    
    i = model.parent(j);
    while i > 0
        ii = model.vinds{i};
        
        C(ii,jj) = S{i}.' * f1;
        C(jj,ii) = ( Sdot{i}.'*f2 + S{i}.'*f3 ).';
        
        H(ii,jj) = S{i}.' * f2;
        H(jj,ii) = H(ii,jj).';
        
        Hdot(ii,jj) = Sdot{i}.'*f2 + S{i}.'*(f1 + f3);
        Hdot(jj,ii) = Hdot(ii,jj).';
        
        f1 = Xup{i}.' * f1;
        f2 = Xup{i}.' * f2;
        f3 = Xup{i}.' * f3;
        i = model.parent(i);
    end    
    if model.parent(j) ~= 0
        p = model.parent(j);
        
        IC{p} = IC{p} + Xup{j}.'*IC{j}*Xup{j};
        BC{p} = BC{p} + Xup{j}.'*BC{j}*Xup{j};
    end 
end