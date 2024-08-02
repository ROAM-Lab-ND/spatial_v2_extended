function q = normalizeConfVec(model,q)

    if ~isfield(model,'nq')
        model = postProcessModel(model);
    end
    if ~iscell(q)
        [q] = confVecToCell(model,q);
    end
    
    for i = 1:model.NB
        switch model.jtype{i}
            case {'Fb','S'}
                q{i}(1:4) = q{i}(1:4)/norm( q{i}(1:4) );
            otherwise
                q{i} = q{i};
        end
    end
    q = cell2mat(q);
end

        
