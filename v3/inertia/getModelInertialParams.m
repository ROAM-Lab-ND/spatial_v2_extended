function [a, a_rotor] = getModelInertialParams(model)
    if ~isfield(model,'nq')
        model = postProcessModel(model);
    end
    a = [];
    a_rotor = [];
    for i = 1:model.NB
       a = [a ; inertiaMatToVec( model.I{i} )]; 
       if model.has_rotor(i)
          a_rotor = [a_rotor ; inertiaMatToVec( model.I_rotor{i} ) ]; 
       end
    end
end