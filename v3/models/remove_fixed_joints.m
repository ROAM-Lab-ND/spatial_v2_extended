function model = remove_fixed_joints(model)
% Update inertias and transforms
for i = 1:model.NB
    if ~strcmp(model.jtype{i}, 'fixed')
        continue;
    end
    
    parent = model.parent(i);
    children = find(model.parent == i);
    
    % Update the transforms of the children
    for j = children
        model.parent(j) = parent;
        model.Xtree{j} = model.Xtree{j}*model.Xtree{i};
    end
    
    % Absorb the body into the parent body if it exists
    if parent > 0
        % Add inertia to parent
        X = model.Xtree{i};
        I = X'*model.I{i}*X;
        model.I{parent} = model.I{parent} + I;
    end
end

% Remove the joints
i = 1;
while i <= model.NB
    if ~strcmp(model.jtype{i}, 'fixed')
        i = i + 1;
    else
        model.NB = model.NB - 1;
        model.parent(i) = [];
        model.Xtree(i) = [];
        model.I(i) = [];
        model.jtype(i) = [];
        
        for j = i:model.NB
            if model.parent(j) >= i
                model.parent(j) = model.parent(j) - 1;
            end
        end
    end
end
end