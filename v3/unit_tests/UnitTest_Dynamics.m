% Unit tests for dynamics functions

N = 6;

% Create a random model with N links
model = autoTree(N, 1.5, pi/3);
model = postProcessModel(model);
checkDynamics(model,'Fixed Base No Rotors');

model.jtype{1} = 'Fb';
model = postProcessModel(model);
checkDynamics(model,'Floating Base No Rotors');

model.jtype{3} = 'S';
model = postProcessModel(model);
checkDynamics(model,'Floating Base Spherical No Rotors');



model = autoTree_rotor(N, 1.5, pi/3);
model = postProcessModel(model);
checkDynamics(model,'Fixed Base w/ Rotors');

model.jtype{1} = 'Fb';
model.has_rotor(1) = 0;
model = postProcessModel(model);
checkDynamics(model,'Floating Base w/ Rotors');

model.jtype{3} = 'S';
model.has_rotor(3) = 0;
model = postProcessModel(model);
checkDynamics(model,'Floating Base Spherical w/ Rotors');




function checkDynamics(model, desc)
    fprintf('====================================\n');
    fprintf('%s\n',desc);
    fprintf('====================================\n');
    
    % Random inertial properties
    for i = 1:model.NB
        model.I{i} = inertiaVecToMat( rand(10,1) );
        model.I_rotor{i} = inertiaVecToMat( rand(10,1) );
    end
    [a, a_rot] = getModelInertialParams(model);
    
    % Random configuration and velocity
    q   = rand(model.NQ,1);
    q   = normalizeConfVec(model, q); 

    qd  = rand(model.NV,1);
    qd_r= rand(model.NV,1);
    qdd = rand(model.NV,1);

    % Calculate dynamics quanitites
    tau        = ID(model, q ,qd ,qdd);                    % Inverse dynamics
    qdd_ABA    = FDab( model, q, qd, tau);                 % Forward dyanmics
    [H, Cqd]   = HandC(model, q, qd);                      % Mass matrix and bias
    [~, tau_g] = HandC(model, q, qd*0);                    % Gravity force
    Cqd        = Cqd-tau_g;                                % Coriolis force
    [C,Hdot,H2]= CoriolisMatrix( model, q, qd);            % Coriolis matrix
    [Hqd, CTqd]= Hqd_and_CTqd( model, q , qd);             % Gen momentum and friend
    tau_SL     = ID_SlotineLi( model, q, qd , qd_r, qdd);  % Slotine Li ID

    checkValue('H'     , H      , H2                    ); % Mass matrix output is correct
    checkValue('Hqd'   , H*qd   , Hqd                   ); % Generelized momentum
    checkValue('CTqd'  , C'*qd  , CTqd                  ); % For generalized momentum obs
    checkValue('Cqd'   , C*qd   , Cqd                   ); % Generalized Coriolis force
    checkValue('Hdot'  , Hdot   , C+C'                  ); % Hdot -2C skew symmetric
    checkValue('qdd'   , qdd    , qdd_ABA               ); % Forward dynamics
    checkValue('tau'   , tau    , H*qdd+C*qd+tau_g      ); % Inverse dynamics
    checkValue('tau_SL', tau_SL , H*qdd + C*qd_r + tau_g );% Slotine Li inv dyn

    % Make sure that Hdot is correct (finite difference approximation)
    dt = sqrt(eps);
    q_new = q+dt*configurationRates(model,q,qd);
    q_new = normalizeConfVec(model, q_new);
    [H_new, ~] = HandC(model, q_new, qd);
    Hdot_finite_difference = (H_new-H)/dt;
    checkValue('Hdot', Hdot , Hdot_finite_difference, 1e-4 ); % Hdot
    
  

    if ~any(model.has_rotor)
        Hinv = Hinverse(model, q);
        checkValue('Hinv', Hinv , inv(H) );
    end

    % Check Christoffel
    if ~any(model.nv > 1) && ~any(model.has_rotor)
        Gamma = Christoffel(model,q);
        C2 = 0*C;
        for i = 1:model.NB
            C2 = C2 + Gamma(:,:,i)*qd(i);
        end
        Hpartial = H_diff(model,q);
        Gamma2 = 0*Gamma;
        Hdot2 = 0*Hdot;
        for i = 1:model.NB
            for j = 1:model.NB
                for k = 1:model.NB
                    Gamma2(i,j,k) = 1/2* (Hpartial(i,j,k) + Hpartial(i,k,j) - Hpartial(j,k,i));
                end
            end
            Hdot2 = Hdot2 + Hpartial(:,:,i)*qd(i);
        end
        
        checkValue('CGamma'  , C      , C2    ); % Christoffel
        checkValue('Gamma'   , Gamma2 , Gamma ); % Christoffel
        checkValue('Hdot'    , Hdot   , Hdot2 ); % Christoffel
    end
    
    % Regressors
    [Y_Hqd, Y_CTqd, Y_Hqd_rot, Y_CTqd_rot]  ...
                              = Regressor_HqdandCTqd( model, q , qd);
    [Y, Y_rot]                = RegressorClassical(model, q, qd, qdd);
    [Y_SL, Y_SL_rot]          = RegressorSL(model, q, qd,qd_r, qdd);
    [YH, Yg, YH_rot, Yg_rot]  = RegressorHandG(model, q);
    
    if any(model.has_rotor)
        checkValue('Y'   , Y*a    + Y_rot*a_rot    , tau    ); % Classical Regressor
        checkValue('Y_SL', Y_SL*a + Y_SL_rot*a_rot , tau_SL ); % Slotine Li Regressor
        checkValue('YH'  , YH*a   + YH_rot*a_rot   , H(:)   ); % Mass matrix regressor
        checkValue('Yg'  , Yg*a   + Yg_rot*a_rot   , tau_g  ); % Gravity force regressor
        checkValue('YHqd', Y_Hqd*a + Y_Hqd_rot*a_rot , H*qd ); % Indirect regressor
        checkValue('YCTqd',Y_CTqd*a + Y_CTqd_rot*a_rot, C'*qd); % Indirect regressor
    else
        checkValue('Y'     , Y*a      , tau    ); % Classical Regressor
        checkValue('Y_SL'  , Y_SL*a   , tau_SL ); % Slotine Li Regressor
        checkValue('YH'    , YH*a     , H(:)   ); % Mass matrix regressor
        checkValue('Yg'    , Yg*a     , tau_g  ); % Gravity force regressor
        checkValue('Y_Hqd' , Y_Hqd*a  , H*qd   ); % Indirect regressor
        checkValue('Y_CTqd', Y_CTqd*a , C'*qd  ); % Indirect regressor
    end
    
    fprintf('\n');
end

function checkValue(name, v1, v2, tolerance)
    if nargin == 3
        tolerance = sqrt(eps);
    end
    value = norm(v1(:)-v2(:));
    fprintf('%s \t %e\n',name,value);
%     if value > tolerance
%         error('%s is out of tolerance',name);
%     end
end