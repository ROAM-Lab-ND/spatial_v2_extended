
% Unit tests for derivatives functions

clear
N = 3;

% Create a random model with N links
model = autoTree(N, 1.5, pi/3);
checkDerivatives(model,'Floating Base No Rotors');


function checkDerivatives(model, desc)
    fprintf('====================================\n');
    fprintf('%s\n',desc);
    fprintf('====================================\n');
    
    model.jtype{1} = 'Fb';
    model = postProcessModel(model);
    
    % Random inertial properties
    for i = 1:model.NB
        model.I{i} = inertiaVecToMat( rand(10,1) );
    end
    [a] = getModelInertialParams(model);
    
    % Random configuration and velocity
    q   = rand(model.NQ,1);
    if strcmp(model.jtype{1},'SE3')
        q(1:16) = reshape( randomSE3(),[16 1]);
    elseif strcmp(model.jtype{1},'Fb')
        q(1:4) = q(1:4)/norm(q(1:4));
    end
    
    %q   = normalizeConfVec(model, q); 

    qd  = rand(model.NV,1);
    qdd = rand(model.NV,1);
    lambda = rand(model.NV,1);

    % Calculate dynamics quanitites
    out        = modID(model,q,qd,qdd,lambda);
    [tau]      = ID(model, q ,qd ,qdd);                    % Inverse dynamics
    qdd_ABA    = FDab( model, q, qd, tau);                 % Forward dyanmics
    [H, Cqd]   = HandC(model, q, qd);                      % Mass matrix and bias
    [~, tau_g] = HandC(model, q, qd*0);                    % Gravity force
    Cqd        = Cqd-tau_g;                                % Coriolis force
    [C,Hdot,H2]= CoriolisMatrix( model, q, qd);            % Coriolis matrix
    
    checkValue('Cqd'   , C*qd   , Cqd                   ); % Generalized Coriolis force
    checkValue('Hdot'  , Hdot   , C+C'                  ); % Hdot -2C skew symmetric

    % Check Christoffel
    newConfig = @(x) configurationAddition(model,q,x);
    
    dtau_dqd_cs = complexStepJacobian(@(x) ID(model, q ,x ,qdd), qd);
    dtau_dq_cs = complexStepJacobian(@(x) ID(model, newConfig(x) ,qd ,qdd), zeros(model.NV,1) );
    
    dtau_dqd_fd = finiteDiffJacobian(@(x) ID(model, q ,x ,qdd), qd);
    dtau_dq_fd = finiteDiffJacobian(@(x) ID(model, newConfig(x) ,qd ,qdd), zeros(model.NV,1) );
    
    dmodID_dq_cs  = complexStepJacobian(@(x) modID(model, newConfig(x) ,qd ,qdd,lambda), zeros(model.NV,1) );
    dmodID_dqd_cs = complexStepJacobian(@(x) modID(model, q ,x  ,qdd,lambda), qd);
% 
    dmodFD_dq_cs  = complexStepJacobian(@(x) modFD(model, newConfig(x) ,qd ,tau,lambda), 0*qd );
    dmodFD_dqd_cs = complexStepJacobian(@(x) modFD(model, q ,x  ,tau,lambda), qd);
    dmodFD_dtau_cs = complexStepJacobian(@(x) modFD(model,q ,qd ,x ,lambda), tau);

    dqdd_dq_cs  = complexStepJacobian( @(x) FDab(model,newConfig(x) ,qd,tau),zeros(model.NV,1));
    dqdd_dqd_cs = complexStepJacobian( @(x) FDab(model,q,x,tau),qd );
    dqdd_dtau_cs= complexStepJacobian( @(x) FDab(model,q,qd,x),tau );
    
    dqdd_dq_fd  = finiteDiffJacobian( @(x) FDab(model,newConfig(x) ,qd,tau),zeros(model.NV,1));
    dqdd_dqd_fd = finiteDiffJacobian( @(x) FDab(model,q,x,tau),qd );
    dqdd_dtau_fd= finiteDiffJacobian( @(x) FDab(model,q,qd,x),tau );
    

    [dtau_dq, dtau_dqd] = ID_derivatives( model, q, qd, qdd );
    [dqdd_dq, dqdd_dqd,dqdd_dtau] = FD_derivatives( model, q, qd, tau );
    
    [dmodID_dq, dmodID_dqd] = modID_derivatives( model, q, qd, qdd, lambda );
    [dmodFD_dq, dmodFD_dqd, dmodFD_dtau] = modFD_derivatives( model, q, qd, tau, lambda );
% 
    modID_qq_cs   = complexStepJacobian( @(x) outputSelect(1,@modID_derivatives,model,newConfig(x),qd,qdd,lambda), zeros(model.NV,1) );
    modID_qdqd_cs = complexStepJacobian( @(x) outputSelect(2,@modID_derivatives,model,q,x,qdd,lambda),qd );
    modID_qdq_cs  = complexStepJacobian( @(x) outputSelect(2,@modID_derivatives,model,newConfig(x),qd,qdd,lambda), zeros(model.NV,1) );
% 
    derivs = modID_second_derivatives( model, q, qd, qdd, lambda);
    modID_qq = derivs.dmod_dqq;
    modID_qdqd = derivs.dmod_dvv;
    modID_qdq  = derivs.dmod_dqv';
    
    dmod_dq_mid2nd  = derivs.dmod_dq;
    dmod_dqd_mid2nd = derivs.dmod_dv;
    
    dtau_dq_mid2nd  = derivs.dtau_dq;
    dtau_dqd_mid2nd = derivs.dtau_dv;
    
    
%     Hbig  = multiComplexStepHessian(@(x) modID(model, newConfig(x(1:model.NV)) ,x(model.NV+1:end) ,qdd,lambda), [zeros(model.NV,1) ;qd] );

    [~, ~, ~, modFD_qq, modFD_qdqd, modFD_qdq, modFD_tauq] = modFD_second_derivatives( model, q, qd, tau, lambda );
% 
    modFD_qq_cs   = complexStepJacobian( @(x) outputSelect(1,@modFD_derivatives,model,newConfig(x),qd,tau,lambda),0*qd );
    modFD_qdqd_cs = complexStepJacobian( @(x) outputSelect(2,@modFD_derivatives,model,q,x,tau,lambda),qd );
    modFD_qdq_cs  = complexStepJacobian( @(x) outputSelect(2,@modFD_derivatives,model,newConfig(x),qd,tau,lambda),0*qd );
    modFD_tauq_cs = complexStepJacobian( @(x) outputSelect(3,@modFD_derivatives,model,newConfig(x),qd,tau,lambda),0*qd );

    checkValue('ID_q'   , dtau_dq      , dtau_dq_cs   ); % Partials of ID w.r.t. q
    checkValue('ID_qd'  , dtau_dqd     , dtau_dqd_cs  ); % Partials of ID w.r.t. qd
% 
    checkValue('FD_q'   , dqdd_dq      , dqdd_dq_cs   ); % Partials of FD w.r.t. q
    checkValue('FD_qd'  , dqdd_dqd     , dqdd_dqd_cs  ); % Partials of FD w.r.t. qd
    checkValue('FD_tau' , dqdd_dtau    , dqdd_dtau_cs ); % Partials of FD w.r.t. tau
    
    fprintf('Finite Diff ====================================\n');
    
    checkValue('ID_q_fd'   , dtau_dq_fd      , dtau_dq_cs     , 5e-03       ); % Partials of ID w.r.t. q
    checkValue('ID_qd_fd'  , dtau_dqd_fd     , dtau_dqd_cs    , 5e-03       ); % Partials of ID w.r.t. qd

    checkValue('FD_q_fd'   , dqdd_dq_fd      , dqdd_dq_cs     , 5e-03       ); % Partials of FD w.r.t. q
    checkValue('FD_qd_fd'  , dqdd_dqd_fd     , dqdd_dqd_cs    , 5e-03       ); % Partials of FD w.r.t. qd
    checkValue('FD_tau_fd'   , dqdd_dtau_fd    , dqdd_dtau_cs , 5e-03       ); % Partials of FD w.r.t. tau
% 
    disp('Complex Step ====================================');
%   
    checkValue('modID_q'   , dmodID_dq'      , dmodID_dq_cs            ); % Partials of modID w.r.t. q
    checkValue('modID_qd'  , dmodID_dqd'     , dmodID_dqd_cs           ); % Partials of modID w.r.t. qd
% 
    checkValue('modFD_q'   , dmodFD_dq      , dmodFD_dq_cs            ); % Partials of modFD w.r.t. q
    checkValue('modFD_qd'  , dmodFD_dqd     , dmodFD_dqd_cs           ); % Partials of modFD w.r.t. qd
    checkValue('modFD_tau' , dmodFD_dtau    , dmodFD_dtau_cs          ); % Partials of modFD w.r.t. qd
% 
    checkValue('modID_qq'    , modID_qq      , modID_qq_cs            ); % SO Partials of modID w.r.t. q,q
    checkValue('modID_qdq'   , modID_qdq     , modID_qdq_cs           ); % SO Partials of modID w.r.t. qd,q
    checkValue('modID_qdqd'  , modID_qdqd    , modID_qdqd_cs          ); % SO Partials of modID w.r.t. qd,qd
    
    
    checkValue('mod_q'   , dmod_dq_mid2nd    , dmodID_dq_cs           ); % SO Partials of modID w.r.t. qd,qd
    checkValue('mod_qd'  , dmod_dqd_mid2nd   , dmodID_dqd_cs          ); % SO Partials of modID w.r.t. qd,qd
    
    checkValue('ID_q'   , dtau_dq_mid2nd    , dtau_dq_cs           ); % SO Partials of modID w.r.t. qd,qd
    checkValue('ID_qd'  , dtau_dqd_mid2nd   , dtau_dqd_cs          ); % SO Partials of modID w.r.t. qd,qd
    
% 
    checkValue('modFD_qq'    , modFD_qq      , modFD_qq_cs            ); % SO Partials of modID w.r.t. q,q
    checkValue('modFD_qdq'   , modFD_qdq     , modFD_qdq_cs           ); % SO Partials of modID w.r.t. qd,q
    checkValue('modFD_qdqd'  , modFD_qdqd    , modFD_qdqd_cs          ); % SO Partials of modID w.r.t. qd,qd
    checkValue('modFD_tauq'  , modFD_tauq    , modFD_tauq_cs           ); % SO Partials of modID w.r.t. tau,q

    fprintf('\n');
end



function checkValue(name, v1, v2, tolerance)
    if nargin == 3
        tolerance = sqrt(eps);
    end
    value = norm(v1(:)-v2(:));
    fprintf('%10s \t %e\n',name,value);
    if value > tolerance
        error('%s is out of tolerance',name);
    end
end