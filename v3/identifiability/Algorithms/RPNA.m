function [N, M, V , C] = RPNA(model, free_base)
%% [N, M, V , C] = RPNA(model, free_base)
% Recursive parameter nullspace algorithm
% The outputs of the algorithm are
% N: cell array for the nullspace descriptors characterizing undetectable
%    inertial transfers across each joint
% M: Basis for minimal paramers (sometimes called representative base
%    parameters) of each link
% V: Attainable velocity span basis for each link
% C: Null(C) gives the unidentifiable parameters of each link
%
% See the accompanying paper for further detail.
%
    if nargin == 1
        free_base = 0;
    end

    % Initialize Empty Cell arrays
    EmptyCells = cell(model.NB,1); 
    N   = EmptyCells;   C   = EmptyCells;   O   = EmptyCells;
    M   = EmptyCells;   V   = EmptyCells;   
    V_J = EmptyCells;  
    
    % Main algorithm loop
    for i =1:model.NB
        [~, Si] = jcalc( model.jtype{i}, 0 );
        n_i = size(Si, 2);
        
        p = model.parent(i);
        if p == 0
            % If the base is to be treated as free, set the motion and
            % outer product bases to represent the full space. 
            % Else initialize them from the fixed base.
            if free_base
                Vp = eye(6);
                Cp = eye(10);
            else
                Vp =  get_gravity(model);
                Cp = zeros(0,10);
            end
        else
            Vp = V{p};
            Cp = C{p};
        end
        
        % Propage the velocity and outer product spans across the link
        V_J{i} = model.Xtree{i} * Vp;
       
        % Compute the collection of rate matricies for multi-DoF joints
        crmSet = cell( n_i, 1);
        paramRateSet     = cell( n_i, 1);
        for k= 1:size(Si,2)
            % Velocity rates of change with joint angle
            crmSet{k} = crm(Si(:,k) );
            
            % Inertail parameters rates of change
            paramRateSet{k} = Rate_Parameters( Si(:,k) );
        end
        
        % Propagate the velocity across a joint
        V{i} = RangeBasis([ SwitchedControllabilityMatrix( crmSet , V_J{i} ) Si]); 
        
        O{i} = SwitchedObservabilityMatrix(Cp*Transform_Parameters( model.Xtree{i} ), paramRateSet );
        N{i} = OutputMatrix(V{i},Si);
        for k = 1:n_i
            N{i} = [N{i} ; O{i} * paramRateSet{k}];
        end
        N{i} = RangeBasis(N{i}')';
        C{i} = RangeBasis( [ Cp*Transform_Parameters( model.Xtree{i} ) ; N{i} ]')'; 
        M{i} = UnitVectorComplementarySubspace( null(N{i}) );
    end
end