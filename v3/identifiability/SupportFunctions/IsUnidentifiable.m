function bool = IsUnidentifiable(model,N, i, k)
%% bool = IsUnidentifiable(model,N, i, k)
%  Determine if parameter k of body i is unidentifiable
%  Requires the nullspace descriptor array N from the RPNA

    pi = zeros(10,1);
    pi(k) = 1;
    bool = true;
    
    while i > 0
        if norm( N{i}* pi) > eps^.75
            bool = false;
            return
        end
        X  = model.Xtree{i} ;
        pi = inertiaMatToVec( X' * inertiaVecToMat(pi) * X ) ;
        i = model.parent(i);
    end
    