function bool = IsIdentifiable_motor(model,N, i, k)
%% bool = IsIdentifiable_motor(model,N, i, k)
% Determine if parameter k of body i is identifiable
%  if i > model.NB, then body i refers to motor (i - model.NB) 
% Requires the nullspace descriptor array N from the RPNA

    if i <= model.NB
        R = null(N{i});
        pi = zeros(20,1); 
        pi(k) = 1;

        if norm( pi' * R ) > eps^.75
            bool = false;
            return
        end

        for j = 1:model.NB
            pi = pi(1:10);
            if i == model.parent(j)
                Rj = null(N{j});

                AX = Transform_Parameters( model.Xtree{j} );
                AX_rotor = Transform_Parameters( model.Xrotor{j} );

                if norm( pi'*([AX AX_rotor*model.rotor_constraint{i}]*Rj) ) > eps^.75
                    bool=false;
                    return
                end
            end
        end

        bool = true;
    else
        R = null(N{i-model.NB});
        pi = zeros(20,1); 
        pi(10+k) = 1;

        if norm( pi' * R ) > eps^.75
            bool = false;
            return
        end
        bool = true;
    end
end
