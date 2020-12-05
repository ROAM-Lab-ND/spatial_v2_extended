clear all; clc;


%% Options
% Uncomment to pick a model
model = Puma560_model();      % Industrial robot, 6 DoF
% model = Scara_model();        % Industrial robot, 3 DoF
% model = CheetahLeg_model();   % Single leg of a quadruped, 3 DoF

MODEL_MOTORS = 0; % Include motor inertias (1), or ignore them (0)
FIXED_BASE   = 1; % Treat as fixed base (1), or ignore motion restrictions (0)

%model.gravity = [0 0 0]';
num_regressor_samples = 20; % For comparisson to numerical SVD

%% Compute Parameter Nullspace with SVD
fprintf(1,'\n\n***********************************************\n\n');
fprintf('Computing Random Regressors\n')
if MODEL_MOTORS
    Ystack = ComputeSampledRegressor(model, num_regressor_samples);
    for i = 1:model.NB
        sub = Ystack(:,10*(model.NB+i-1)+1: 10*(model.NB+i-1)+10) ;
        Ystack(:,10*(model.NB+i-1)+1: 10*(model.NB+i-1)+10) = sub*model.rotor_constraint{i};
    end
else
    model.has_rotor = zeros(model.NB,1);
    Ystack = ComputeSampledRegressor(model, num_regressor_samples);
end
[Uy, Ey, Vy] = svd(Ystack);

if MODEL_MOTORS
    params_per_body = 20;
else
    params_per_body = 10;
end
SVD_Nullspace_Dimension = params_per_body*model.NB - rank(Ystack);

%RangeBasis(1,1);
%[~,SVD_Condition] = RangeBasis(Ystack');

%% Compute Parameter Nullspace with RPNA
% RangeBasis(1,1);
param_names = {'m', 'mcx', 'mcy', 'mcz', 'Ixx', 'Iyy', 'Izz', 'Iyz', 'Ixz', 'Ixy'};
fprintf('Running RPNA\n');
fprintf(1,'\n===============================\n');
fprintf(1,'Identifiable Parameter Detail\n');

if MODEL_MOTORS
    [N, M, V, C] = RPNA_rotor(model,~FIXED_BASE);
    %[~, RPNA_Condition] = RangeBasis(1);
    [Null_Basis, Minimal_Basis, Perp_Basis, Perp_Basis_sym] = ComputeBases_rotor(model, N, M);
    PrintParameterSummary_rotor(model, N, M, V,C, param_names);
else
    [N, M, V, C] = RPNA(model,~FIXED_BASE);
    %[~, RPNA_Condition] = RangeBasis(1);
    [Null_Basis, Minimal_Basis, Perp_Basis, Perp_Basis_sym] = ComputeBases(model, N, M);
    PrintParameterSummary(model, N, M, V,C, param_names);
end

RPNA_Nullspace_Dimension = 0;
for i = 1:model.NB
    RPNA_Nullspace_Dimension = RPNA_Nullspace_Dimension + params_per_body-size(N{i},1);
end


%% Compute identifiable linear combinations with rref
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



% Create variables for printing parameter regroupings
sym_params = sym( zeros(params_per_body*model.NB,1) ) ;    
for i = 1:model.NB
    for k = 1:10
        sym_params(10*i-10 + k ) = sym(sprintf('%s%d',param_names{k},i));
        if MODEL_MOTORS
            sym_params(10*(i+model.NB)-10 + k ) = sym(sprintf('%sM%d',param_names{k},i));
        end
    end
end

% Compute identifiable parameter combinations from the basis for the
% subspace perpendicular to the parameter nullspace
Perp_Basis = rref(Perp_Basis')';
Perp_Basis_sym = rref(Perp_Basis_sym')';

inds = find(abs(Perp_Basis) < 1e-8); % remove small values so printing is clean
Perp_Basis(inds) = 0;
Perp_Basis_sym(inds) = 0;

inds = find(abs(Perp_Basis-1) < 1e-8); % remove small values so printing is clean
Perp_Basis(inds) = 1;
Perp_Basis_sym(inds) = 1;

inds = find(abs(Perp_Basis+1) < 1e-8); % remove small values so printing is clean
Perp_Basis(inds) = -1;
Perp_Basis_sym(inds) = -1;

regrouping_matrix = sym(zeros(params_per_body*model.NB, params_per_body*model.NB  ));
for i = 1:size(Perp_Basis_sym,2)
    ind = find(Perp_Basis_sym(:,i)==1,1);
    regrouping_matrix(ind, :) = Perp_Basis_sym(:,i)';
    
    % Identifable parameter combination
    sym_result = Perp_Basis_sym(:,i)'*sym_params;
    
    % Work to strip out zero coefficients
    [coef, monomials] = coeffs(sym_result);
    coef = CleanMat(coef);
    sym_result = simplify( sum( coef(:).*monomials(:) ) );
    
    % And then group terms that multiply each parameter
    sym_result = jacobian(sym_result, sym_params)*sym_params;
    fprintf(1,'Regrouped parameter %s <= %s\n', char(sym_params(ind)), char(sym_result));
end

fprintf(1,'\n===================================\n');
fprintf(1,'Sanity Checks \n');
fprintf(1,'===================================\n');
fprintf('Null Check => norm( Ystack * Null_Basis ) = %e\n', norm( Ystack * Null_Basis , 'fro'))
fprintf('Perp Check => norm( Null_Basis\''*Perp_Basis ) = %e \n\n',norm(Null_Basis'*Perp_Basis,'fro') );

fprintf(1,'===================================\n');
fprintf(1,'Summary \n');
fprintf(1,'===================================\n');

if FIXED_BASE
    fprintf('Nullspace Dimension SVD  %d\n',SVD_Nullspace_Dimension)
end
fprintf('Nullspace Dimension RPNA %d\n',RPNA_Nullspace_Dimension)
fprintf('Identifiable Dimension %d\n',model.NB*params_per_body - RPNA_Nullspace_Dimension)

fprintf(1,'\n');