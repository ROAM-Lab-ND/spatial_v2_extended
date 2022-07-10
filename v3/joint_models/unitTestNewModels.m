clear
%% Test with rotors
disp('Absolute Pair Joints with Rotors')
org_model = autoTree_rotor(8,2,.3);

Xtree = org_model.Xtree;

model =  RBD_model();
jointAxes = {'x','y','z','x','y','z','x','y','z','x','y','z'};

model.NB = org_model.NB;
model.parent = org_model.parent;


model.Xtree = org_model.Xtree;
model.I     = org_model.I;
    

model.joint{1} = floatingBaseJoint();
model.I{1} = randomInertia();

a = [];
a(1:10,1) = inertiaMatToVec( model.I{1} );
for i = 2:model.NB
    model.joint{i} = revolutePairAbsoluteWithRotors();
    model.joint{i}.jointAxis{1} = jointAxes{i};  % x, y, or z
    model.joint{i}.jointAxis{2} = jointAxes{i};  % x, y, or z
    
    model.joint{i}.rotorAxis{1} = jointAxes{i}; % x, y, or z
    model.joint{i}.rotorAxis{2} = jointAxes{i}; % x, y, or z
    
    model.joint{i}.XtreeInternal = randomAdSE3(); % from first link to second link
    
    model.Xtree{i}(1:6,:) = Xtree{i}; % from predecessor to first link
    model.Xtree{i}(7:12,:) = org_model.Xrotor{i}; % from predecessor to first rotor
    model.Xtree{i}(13:18,:) = randomAdSE3(); % from predecessor to second rotor
    
    model.I{i}(1:6,1:6) = randomInertia(); % first link (shank)
    model.I{i}(7:12,7:12) = randomInertia(); % first rotor
    model.I{i}(13:18,13:18) = randomInertia(); % second rotor
    model.I{i}(19:24,19:24) = randomInertia(); % second link (foot)
    
    
    a(end+1 : end+10) = inertiaMatToVec(model.I{i}(1:6,1:6));
    a(end+1 : end+10) = inertiaMatToVec(model.I{i}(7:12,7:12));
    a(end+1 : end+10) = inertiaMatToVec(model.I{i}(13:18,13:18));
    a(end+1 : end+10) = inertiaMatToVec(model.I{i}(19:24,19:24));
    
    model.joint{i}.gearRatio{1} = rand()*5; % first rotor
    model.joint{i}.gearRatio{2} = rand()*5; % second rotor
    model.joint{i}.beltRatio{1} = rand()*5;
    model.joint{i}.beltRatio{2} = rand()*5;

%% Rotors
%     model.joint{i} = revoluteJointWithRotor();
%     model.joint{i}.jointAxis = jointAxes{i};
%     model.joint{i}.rotorAxis = jointAxes{i};
%     model.joint{i}.gearRatio = rand()*5;
%     model.I{i}(1:6,1:6) = randomInertia(); % first link (shank)
%     model.I{i}(7:12,7:12) = randomInertia(); % first rotor
%     
%     model.Xtree{i}(1:6,:) = Xtree{i}; % from predecessor to first link
%     model.Xtree{i}(7:12,:) = org_model.Xrotor{i};
% 
% 
%     a(end+1 : end+10) = inertiaMatToVec(model.I{i}(1:6,1:6));
%     a(end+1 : end+10) = inertiaMatToVec(model.I{i}(7:12,7:12));
end

model = model.postProcessModel();

q = rand(model.NQ,1); % using q's for the joints as relative link angles

q = model.normalizeConfVec(q);
qd = rand(model.NV,1);
qdr = rand(model.NV,1);
qdd = rand(model.NV,1);

% Random configuration and velocity
q   = rand(model.NQ,1);
q   = normalizeConfVec(model, q); 

qd  = rand(model.NV,1);
qd_r= rand(model.NV,1);
qdd = rand(model.NV,1);
lambda=rand(model.NV,1);

% Calculate dynamics quanitites
[tau, out] = ID(model, q ,qd ,qdd);                    % Inverse dynamics
qdd_ABA    = FDab( model, q, qd, tau);                 % Forward dyanmics
[H, Cqd, info] = HandC(model, q, qd);                      % Mass matrix and bias
p          = H*qd;
[~, tau_g] = HandC(model, q, qd*0);                    % Gravity force
Cqd        = Cqd-tau_g;                                % Coriolis force
[C,Hdot,H2]= CoriolisMatrix( model, q, qd);            % Coriolis matrix
[Hqd, CTqd]= Hqd_and_CTqd( model, q , qd);             % Gen momentum and friend
tau_SL     = ID_SlotineLi( model, q, qd , qd_r, qdd);  % Slotine Li ID
[qdot, pdot] = HamiltonianDynamics(model,q,p,tau);

checkValue('H'     , H      , H2                    ); % Mass matrix output is correct
checkValue('Hqd'   , H*qd   , Hqd                   ); % Generelized momentum
checkValue('CTqd'  , C'*qd  , CTqd                  ); % For generalized momentum obs
checkValue('Cqd'   , C*qd   , Cqd                   ); % Generalized Coriolis force
checkValue('Hdot'  , Hdot   , C+C'                  ); % Hdot -2C skew symmetric
checkValue('qdd'   , qdd    , qdd_ABA               ); % Forward dynamics
checkValue('tau'   , tau    , H*qdd+C*qd+tau_g      ); % Inverse dynamics
checkValue('tau_SL', tau_SL , H*qdd + C*qd_r + tau_g );% Slotine Li inv dyn

checkValue('Ham_qd', qdot , qd );% Slotine Li inv dyn

ret = EnerMo(model, q, qd);
checkValue('Kin'     , ret.KE      , 1/2*qd'*H*qd   ); % Mass matrix output is correct

if strcmp(class(model.joint{1}), 'floatingBaseJoint')
    X1 = out.Xup{1};
    checkValue('Itot'     , ret.Itot      , X1'*info.IC{1}*X1   ); % Mass matrix output is correct
    checkValue('htot'     , ret.htot      , X1'*H(1:6,:)*qd);
end


% Make sure that Hdot is correct (complex step approximation)
dt = sqrt(eps)*1i;
q_new = configurationAddition(model,q,dt*qd);
%q_new = normalizeConfVec(model, q_new);
qd_new = qd + dt*qdd;

[H_new, ~] = HandC(model, q_new, qd);
p_new = H_new*qd_new;

Hdot_finite_difference = real( (H_new-H) / dt );
pdot_finite_difference = real( (p_new-p) / dt );

checkValue('Hdot', Hdot , Hdot_finite_difference ); % Hdot
checkValue('pdot', pdot , pdot_finite_difference ); % Hdot

Hinv = Hinverse(model, q);
checkValue('Hinv', Hinv , inv(H) );

out = modID(model,q,qd,qdd,lambda);
checkValue('modID',out, lambda'*tau);

% Check Christoffel

approved_model = 1;
compatible_joints = {'floatingBaseJoint','revoluteJoint','revoluteJointWithRotor'};

for i = 1:model.NB
    if ~any( strcmp(compatible_joints, class(model.joint{i}) ))
        approved_model = 0;
    end
end

    
if approved_model

    Gamma = Christoffel(model,q);

    dC_dqd = complexStepJacobian( @(x) CoriolisMatrix( model, q, x) , 0*qd);
    Gamma_cs = permute(dC_dqd,[1 3 2]);
    checkValue('Gam_cs'  , Gamma_cs      , Gamma ); % Christoffel


    C2 = 0*C;
    for i = 1:model.NV
        C2 = C2 + squeeze(Gamma(:,i,:)*qd(i));
    end
    checkValue('CGamma'  , C      , C2    ); % Christo

    newConfig = @(x) configurationAddition(model,q,x);
    H_partial_cs = complexStepJacobian( @(x) HandC(model, newConfig(x),qd), 0*qd);

    Hpartial = H_derivatives(model,q);
    checkValue('dH_cs'  , Hpartial      , H_partial_cs    ); % Christoffel

    struc_const = StructureConstants(model,q);

    Gamma_via_kozul = 0*Gamma;
    for i = 1:model.NV
        for j = 1:model.NV
            for k = 1:model.NV
                Gamma_via_kozul(i,j,k) = ...
                       Hpartial(k,i,j) + Hpartial(j,i,k) - Hpartial(j,k,i) ...
                       + struc_const(i,j,k) + struc_const(k,i,j) + struc_const(j,i,k);
            end
        end
    end
    Gamma_via_kozul = Gamma_via_kozul/2;

    checkValue('Struc'    , Gamma_via_kozul   , Gamma ); % Christoffel

    Hdot2 = 0*Hdot;
    for i = 1:model.NV
        Hdot2 = Hdot2 + Hpartial(:,:,i)*qd(i);
    end
    checkValue('Hdot'    , Hdot   , Hdot2 ); % Christoffel
end


% Regressors
[Y_Hqd, Y_CTqd] = Regressor_HqdandCTqd( model, q , qd);
[Y]             = RegressorClassical(model, q, qd, qdd);
[Y_SL]          = RegressorSL(model, q, qd,qd_r, qdd);
[YH, Yg]        = RegressorHandG(model, q);


checkValue('Y'     , Y*a      , tau    ); % Classical Regressor
checkValue('Y_SL'  , Y_SL*a   , tau_SL ); % Slotine Li Regressor
checkValue('YH'    , YH*a     , H(:)   ); % Mass matrix regressor
checkValue('Yg'    , Yg*a     , tau_g  ); % Gravity force regressor
checkValue('Y_Hqd' , Y_Hqd*a  , H*qd   ); % Indirect regressor
checkValue('Y_CTqd', Y_CTqd*a , C'*qd  ); % Indirect regressor


J = BodyJacobian( model, q, 8, [zeros(6,18) eye(6)]);
Lambda_inv = J*(H\J');
F_spatial = [1 0 0 0 0 0]';
F_group = [zeros(18,1) ; F_spatial];
[qdd, lambda_inv_entry, sub, D] = applyTestForce(model,q,8,F_group);
qdd2 = H\(J'*F_spatial);

checkValue('Lambda_inv', F_spatial'*Lambda_inv*F_spatial,  lambda_inv_entry)
checkValue('qdd_test', qdd,  qdd2)
checkValue('invHFactor', sub*(D\(sub')), inv(H) ); 



function checkValue(name, v1, v2, tolerance)
    if nargin == 3
        tolerance = sqrt(eps);
    end
    value = norm(v1(:)-v2(:));
    fprintf('%s \t %e\n',name,value);
    if value > tolerance
        error('%s is out of tolerance',name);
    end
end








