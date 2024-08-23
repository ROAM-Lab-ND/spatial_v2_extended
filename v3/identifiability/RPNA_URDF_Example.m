clear all; clc;

% Import model via URDF
floating_base_flag = 0; % Fixed base robot
model = URDF_to_spatialv2_model('panda_arm_no_fixed.urdf');
model.gravity = [0 0 -9.81]';

%% Compute Parameter Nullspace with RPNA
fprintf('Running RPNA\n');
[N, M, V, C] = RPNA(model, floating_base_flag);
[Null_Basis, Minimal_Basis, Perp_Basis, Perp_Basis_sym] = ComputeBases(model, N, M);

% Compute identifiable parameter combinations from the basis for the
% subspace perpendicular to the parameter nullspace
Perp_Basis_sym = rref(Perp_Basis_sym')';
Perp_Basis = rref(Perp_Basis')';

%% Print out parameter regroupings
fprintf(1,'===================================\n');
fprintf(1,'Minimal Parameter Detail\n');
fprintf(1,'===================================\n');
fprintf('Note 1: The listed linear cobminations of parameters are identifiable\n');
fprintf('from fully exciting data. These regroupings are also called minimal\n');
fprintf('parameters or base parameters in the literature. \n\n');

fprintf('Note 2: More geometrically, the regrouping coefficients for each \n')
fprintf('base parameter provide a basis vector for the orthogonal complement to \n')
fprintf('the parameter nullspace \\mathcal{N}. Since the choice of a basis for \n')
fprintf('a vector subspace is not unique, it follows that the choice of base\n'); 
fprintf('parameters is not unique either\n\n');

fprintf('Note 3: The URDF has a separate link frame at each link CoM and with two frames \n')
fprintf('implied around each joint. If Joint i connects body i to its parent, then frame i+ \n')
fprintf('is immediately after the joint and is connected to body i, with i- immediately\n')
fprintf('before the joint and connected to the parent body. All quanities in spatial_v2 \n')
fprintf('are given in the i+ frames. As such the regrouped parameters are given relative \n')
fprintf('to these frames as well.\n\n')

param_names = {'m', 'mcx', 'mcy', 'mcz', 'Ixx', 'Iyy', 'Izz', 'Iyz', 'Ixz', 'Ixy'};

% Create variables for printing parameter regroupings
sym_params = sym( zeros(10*model.NB,1) ) ;    
for i = 1:model.NB
    for k = 1:10
        sym_params(10*i-10 + k ) = sym(sprintf('%s%d',param_names{k},i));
    end
end

sympref('FloatingPointOutput',true);
for i = 1:size(Perp_Basis_sym,2)
    ind = find(Perp_Basis_sym(:,i)==1,1);
    
    % Identifable parameter combination
    sym_result = Perp_Basis_sym(:,i)'*sym_params;
    
    % Work to strip out zero coefficients
    [coef, monomials] = coeffs(sym_result);
    coef = CleanMat(coef);
    sym_result = simplify( sum( coef(:).*monomials(:) ) );
    
    % And then group terms that multiply each parameter
    sym_result = jacobian(sym_result, sym_params)*sym_params;
    fprintf(1,'Regrouped parameter %s <= ', char(sym_params(ind)));
    disp(sym_result)
end

sympref('FloatingPointOutput','default');
fprintf('\nNullspace Dimension RPNA %d\n',size(Null_Basis,2))
fprintf('Identifiable Dimension %d\n',size(Perp_Basis,2))

fprintf(1,'\n');