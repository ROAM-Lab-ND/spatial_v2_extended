function bool = IsIdentifiable(model,N, i, k)
%% bool = IsIdentifiable(model,N, i, k)
% Determine if parameter k of body i is identifiable
% Requires the nullspace descriptor array N from the RPNA

    R = null(N{i});
    pi = zeros(10,1); 
    pi(k) = 1;
    
    if norm( pi' * R ) > eps^.75
        bool = false;
        return
    end
    
    for j = 1:model.NB
        if i == model.parent(j)
            Rj = null(N{j});
            AX = Transform_Parameters( model.Xtree{j} );
            
            if norm( pi'*(AX*Rj) ) > eps^.75
                bool=false;
                return
            end
        end
    end
    
    bool = true;
end