function  [C, Hdot, H ] = CoriolisMatrix( model, q, qd, factorFunction)

if nargin == 3
    factorFunction = @(I,v)(factorFunctions(I,v));
end

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
[q, qd] = confVecToCell(model,q,qd);

I = model.I;
if sum(model.has_rotor) > 0
    I_rotor = model.I_rotor;
end

for i = 1:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
    Xup{i} = XJ * model.Xtree{i};
    if model.parent(i) == 0
        v{i} = S{i}*qd{i};
    else
        v{i} = Xup{i}*v{model.parent(i)} + S{i}*qd{i};
    end
    Sdot{i} = crm(v{i}) * S{i};
    BC{i} = factorFunction(I{i}, v{i});
    
    % Extra data for rotors
    if model.has_rotor(i)
      [ XJ_rotor, S_rotor{i} ] = jcalc( model.jtype_rotor{i}, q{i}*model.gr{i} );
      S_rotor{i} = S_rotor{i} * model.gr{i};
      vJ_rotor = S_rotor{i} * qd{i};
      Xup_rotor{i} = XJ_rotor * model.Xrotor{i};
      if model.parent(i) == 0
          v_rotor{i} = vJ_rotor;
      else
          v_rotor{i} = Xup_rotor{i}*v{model.parent(i)} + vJ_rotor;
      end
      Sdot_rotor{i} = crm(v_rotor{i}) * S_rotor{i};
      
      B_rotor{i} = factorFunction(I_rotor{i}, v_rotor{i});
    end
    
end

IC = I;
C = zeros(model.NV);
H = zeros(model.NV);
Hdot = zeros(model.NV);

for j = model.NB:-1:1
    jj = model.vinds{j};
   
    f1 = IC{j} * Sdot{j} + BC{j} * S{j};
    f2 = IC{j} * S{j};
    f3 = BC{j}' * S{j};
    
    C(jj,jj) = S{j}' * f1;
    H(jj,jj) = S{j}' * f2;
    Hdot(jj,jj) = Sdot{j}'*f2 + S{j}'*(f1 + f3);
    
    f1 = Xup{j}'*f1;
    f2 = Xup{j}'*f2;
    f3 = Xup{j}'*f3;
    
    if model.has_rotor(j)
        f1_rotor = I_rotor{j} * Sdot_rotor{j} + B_rotor{j} * S_rotor{j};
        f2_rotor = I_rotor{j} * S_rotor{j};
        f3_rotor = B_rotor{j}' * S_rotor{j};
        
        C(jj,jj) = C(jj,jj) + S_rotor{j}' * f1_rotor;
        H(jj,jj) = H(jj,jj) + S_rotor{j}' * f2_rotor;
        Hdot(jj,jj) = Hdot(jj,jj) + Sdot_rotor{j}'*f2_rotor + S_rotor{j}'*(f1_rotor + f3_rotor);
    
        f1 = f1 + Xup_rotor{j}'*f1_rotor;
        f2 = f2 + Xup_rotor{j}'*f2_rotor;
        f3 = f3 + Xup_rotor{j}'*f3_rotor;
    end
    
    i = model.parent(j);
    while i > 0
        ii = model.vinds{i};
        
        C(ii,jj) = S{i}' * f1;
        C(jj,ii) = ( Sdot{i}'*f2 + S{i}'*f3 )';
        
        H(ii,jj) = S{i}' * f2;
        H(jj,ii) = H(ii,jj)';
        
        Hdot(ii,jj) = Sdot{i}'*f2 + S{i}'*(f1 + f3);
        Hdot(jj,ii) = Hdot(ii,jj)';
        
        f1 = Xup{i}' * f1;
        f2 = Xup{i}' * f2;
        f3 = Xup{i}' * f3;
        i = model.parent(i);
    end    
    if model.parent(j) ~= 0
        p = model.parent(j);
        
        IC{p} = IC{p} + Xup{j}'*IC{j}*Xup{j};
        BC{p} = BC{p} + Xup{j}'*BC{j}*Xup{j};
        
        if model.has_rotor(j)
            IC{p} = IC{p} + Xup_rotor{j}'*I_rotor{j}*Xup_rotor{j};
            BC{p} = BC{p} + Xup_rotor{j}'*B_rotor{j}*Xup_rotor{j};
        end
    end 
end