function [Null_Basis, Minimal_Basis, Perp_Basis, Perp_Basis_sym] = ComputeBases_rotor(model, N, M) 
%%  [Null_Basis, Minimal_Basis, Perp_Basis, Perp_Basis_sym] = ComputeBases_rotor(model, N, M) 
% Computes system level bases for the parameter nullspace and its
% orthogonal complement. It also gives a minimal basis of unit vectors
% whose span is complementary to the parameter nullspace. Cell array inputs
% N and M are given as outputs of the RPNA algorithm


    num_params = model.NB*20;
    perp_dim = 0; % dimension of the orthogonal complement to the parameter nullspace
    null_dim = 0; % dimension of the parameter nullspace
    perp_inds = cell(model.NB,1); % cell array containing column indicies of the perp_basis for each joint
    null_inds = cell(model.NB,1); % cell array containing column indicies of the null_basis for each joint
    
    for i = 1:model.NB
        dim_i = size(N{i},1);
     
        perp_inds{i} = perp_dim +1 : perp_dim + dim_i;
        null_inds{i} = null_dim +1 : null_dim + 20 - dim_i;
        
        perp_dim = perp_dim + dim_i;
        null_dim = null_dim + 20 - dim_i;
    end
    
    % Matrices whose colums will provide bases for each space
    Null_Basis    = zeros( num_params, null_dim );
    Perp_Basis    = zeros( num_params, perp_dim );
    Perp_Basis_sym= sym(zeros( num_params, perp_dim ));
    Minimal_Basis = zeros( num_params, perp_dim );
    
   
    R    = cell(model.NB,1);  % Holds the inerital transfer basis for each joint
    
    % Loop through one by at a time
    for i = 1:model.NB
        i_inds = parameter_inds(i,model); % Indicies of the inertial parameters for motor & body after joint i
        parent = model.parent(i);   
        
        R{i} =  null(N{i}) ;

        Minimal_Basis( i_inds , perp_inds{i} ) = M{i};
        
        % Undetectable inertial transfers from the child
        Null_Basis   ( i_inds , null_inds{i} ) = R{i};

        if parent > 0
            % X transforms velocities from parent to child
            X = model.Xtree{i};
            X_sym = model.Xtree_sym{i};
            
            % transforms parameters from child to parent
            AX = Transform_Parameters(X);
            AX_sym = Transform_Parameters(X_sym);
            
            AX_rotor = Transform_Parameters(model.Xrotor{i});
            AX_rotor_sym = Transform_Parameters(model.Xrotor_sym{i});
            
            parent_inds = parameter_inds(parent,model); 
            parent_inds = parent_inds(1:10); % only take the indicies for the rigid body of the parent (i.e., don't include its motor)
            
            % This next line gives the undetectable inertial transfers
            % across the joint, as felt by the parent
            Null_Basis(parent_inds, null_inds{i} )= -[AX AX_rotor*model.rotor_constraint{i}]*R{i};
            
            Perp_Basis(i_inds,:)                  = [AX AX_rotor*model.rotor_constraint{i}]'*Perp_Basis(parent_inds,:);
            Perp_Basis_sym(i_inds,:)              = [AX_sym AX_rotor_sym*model.rotor_constraint{i}]'*Perp_Basis_sym(parent_inds,:);
                        
        end
        
        % Set entries close to 0,1,-1 to 0,1,-1 so that the symbolic output
        % will be clean at the end.
        tmp = rref(N{i})';
        inds = find(abs(tmp) < 1e-7); tmp(inds) = 0;      
        inds = find(abs(tmp-1) < 1e-7); tmp(inds) = 1;
        inds = find(abs(tmp+1) < 1e-7); tmp(inds) = -1;
        
        Perp_Basis(i_inds, perp_inds{i} )     = tmp;
        Perp_Basis_sym(i_inds, perp_inds{i} ) = tmp; 
    end
    
    
    %Perp_Basis = rref(Perp_Basis')';
    %Perp_Basis_sym = rref(Perp_Basis_sym')';

    inds = find(abs(Perp_Basis) < 1e-8); % remove small values so printing is clean
    Perp_Basis(inds) = 0;
    Perp_Basis_sym(inds) = 0;

    inds = find(abs(Perp_Basis-1) < 1e-8); % remove small values so printing is clean
    Perp_Basis(inds) = 1;
    Perp_Basis_sym(inds) = 1;

    inds = find(abs(Perp_Basis+1) < 1e-8); % remove small values so printing is clean
    Perp_Basis(inds) = -1;
    Perp_Basis_sym(inds) = -1;

end

function inds = parameter_inds(i,model)
    inds = [(1:10) + 10*(i-1) , (1:10) + 10*(model.NB+i-1)];
end