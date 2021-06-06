function model = postProcessModel(model)
    qi = 0;
    vi = 0;
    ai = 0;
    
    if ~isfield(model,'has_rotor')
        model.has_rotor = zeros(model.NB,1);
    end
    
    
    model.nq = [];
    model.nv = [];
    model.qinds = {};
    model.vinds = {};
    model.ancestors = {};
    model.ancestor_vinds = {};
    model.subtree_vinds = {};
    
    
    for i = 1:model.NB
        [nqi, nvi] = jinfo( model.jtype{i} );
        model.qinds{i} = qi+1 : qi+nqi;
        model.vinds{i} = vi+1 : vi+nvi;
        
        if model.parent(i) == 0
            model.ancestors{i} = [];
            model.ancestor_vinds{i} = [];
        else
            model.ancestors{i} = model.ancestors{ model.parent(i)};
            model.ancestor_vinds{i} = model.ancestor_vinds{ model.parent(i)};
        end
        model.ancestors{i} = [model.ancestors{i}  i];
        model.ancestor_vinds{i} = [model.ancestor_vinds{i} model.vinds{i}];
        
        if model.has_rotor(i)
            model.rotor_param_inds{i} = ai+1  : ai+10;
            ai = ai+10;
        end
        
        model.nq(i) = nqi;
        model.nv(i) = nvi;
        
        qi = qi + nqi;
        vi = vi + nvi;
        model.subtree_vinds{i} = [];
        model.successor_vinds{i} = [];
    end
    
    for i = model.NB:-1:1
        ii = model.vinds{i};
        model.subtree_vinds{i} = [ii model.subtree_vinds{i}];
        if model.parent(i) > 0
            p = model.parent(i);
            model.subtree_vinds{p} = [model.subtree_vinds{i} model.subtree_vinds{p}];
            model.successor_vinds{p} = [ii model.successor_vinds{i} model.successor_vinds{p}];
        end
    end
    
    model.NV = vi;
    model.NQ = qi;
    model.NB_rot = ai/10;
end

