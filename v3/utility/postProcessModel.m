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
    
    for i = 1:model.NB
        [nqi, nvi] = jinfo( model.jtype{i} );
        model.qinds{i} = qi+1 : qi+nqi;
        model.vinds{i} = vi+1 : vi+nvi;
        
        if model.has_rotor(i)
            model.rotor_param_inds{i} = ai+1  : ai+10;
            ai = ai+10;
        end
        
        model.nq(i) = nqi;
        model.nv(i) = nvi;
        
        qi = qi + nqi;
        vi = vi + nvi;
    end
    
    model.NV = vi;
    model.NQ = qi;
    model.NB_rot = ai/10;
end

