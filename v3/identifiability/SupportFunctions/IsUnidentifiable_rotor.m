function bool = IsUnidentifiable_rotor(model,N, i, k)
%% bool = IsUnidentifiable_rotor(model,N, i, k)
% Determine if parameter k of body i is unidentifiable
%  if i > model.NB, then body i refers to motor (i - model.NB) 
% Requires the nullspace descriptor array N from the RPNA

    pi = zeros(20,1);
    if i <= model.NB    
        pi(k) = 1;
    else
        pi(k+10) = 1;
        if norm( N{i-model.NB}* pi) > eps^.75
            bool = false;
            return
        end
        ii = i - model.NB;
        X  = model.Xtree{ii} ;
        X_rotor = model.Xrotor{ii};
        AX = Transform_Parameters( X);
        AX_rotor = Transform_Parameters( X_rotor);
        pi = [ [AX AX_rotor]*pi ; zeros(10,1) ];
        i = i - model.NB;
    end
    
    bool = true;
    while i > 0
        if norm( N{i}* pi) > eps^.75
            bool = false;
            return
        end
        X  = model.Xtree{i} ;
        X_rotor = model.Xrotor{i};
        AX = Transform_Parameters( X);
        AX_rotor = Transform_Parameters( X_rotor);
        pi = [ [AX AX_rotor]*pi ; zeros(10,1) ];
        i = model.parent(i);
    end
    