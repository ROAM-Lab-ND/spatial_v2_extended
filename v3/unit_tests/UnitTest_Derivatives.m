% Unit tests for derivatives functions

clear
N = 15;

% Create a random model with N links
model = autoTree(N, 1.5, pi/3);
model = postProcessModel(model);
checkDerivatives(model,'Fixed Base No Rotors');


function checkDerivatives(model, desc)
    fprintf('====================================\n');
    fprintf('%s\n',desc);
    fprintf('====================================\n');
    
    % Random inertial properties
    for i = 1:model.NB
        model.I{i} = inertiaVecToMat( rand(10,1) );
    end
    [a] = getModelInertialParams(model);
    
    % Random configuration and velocity
    q   = rand(model.NQ,1);
    q   = normalizeConfVec(model, q); 

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
    if ~any(model.nv > 1) && ~any(model.has_rotor)
        Gamma = Christoffel(model,q);
        C2 = 0*C;
        for i = 1:model.NB
            C2 = C2 + Gamma(:,:,i)*qd(i);
        end
        Hpartial = H_derivatives(model,q);
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
        n2 = model.NV*model.NV;
        Hpartial = reshape(Hpartial,n2,model.NV);
        
        dH_dq = complexStepJacobian(@(x) reshape(HandC(model, x, qd),n2,1), q);  
        
        dtau_dqd_cs = complexStepJacobian(@(x) ID(model, q ,x ,qdd), qd);
        dtau_dq_cs = complexStepJacobian(@(x) ID(model, x ,qd ,qdd), q );
        
        dmodID_dq_cs  = complexStepJacobian(@(x) modID(model, x ,qd ,qdd,lambda), q );
        dmodID_dqd_cs = complexStepJacobian(@(x) modID(model, q ,x  ,qdd,lambda), qd);
        
        dmodFD_dq_cs  = complexStepJacobian(@(x) modFD(model, x ,qd ,tau,lambda), q );
        dmodFD_dqd_cs = complexStepJacobian(@(x) modFD(model, q ,x  ,tau,lambda), qd);
        dmodFD_dtau_cs = complexStepJacobian(@(x) modFD(model,q ,qd ,x ,lambda), tau);
        
        
        dqdd_dq_cs  = complexStepJacobian( @(x) FDab(model,x,qd,tau),q );
        dqdd_dqd_cs = complexStepJacobian( @(x) FDab(model,q,x,tau),qd );
        dqdd_dtau_cs= complexStepJacobian( @(x) FDab(model,q,qd,x),tau );
        
        [dtau_dq, dtau_dqd] = ID_derivatives( model, q, qd, qdd );
        [dqdd_dq, dqdd_dqd,dqdd_dtau] = FD_derivatives( model, q, qd, tau );
        [dmodID_dq, dmodID_dqd] = modID_derivatives( model, q, qd, qdd, lambda );
        [dmodFD_dq, dmodFD_dqd, dmodFD_dtau] = modFD_derivatives( model, q, qd, tau, lambda );


        modID_qq_cs   = complexStepJacobian( @(x) outputSelect(1,@modID_derivatives,model,x,qd,qdd,lambda),q );
        modID_qdqd_cs = complexStepJacobian( @(x) outputSelect(2,@modID_derivatives,model,q,x,qdd,lambda),qd );
        modID_qdq_cs  = complexStepJacobian( @(x) outputSelect(2,@modID_derivatives,model,x,qd,qdd,lambda),q );

        derivs = modID_second_derivatives( model, q, qd, qdd, lambda);
        modID_qq = derivs.dmod_dqq;
        modID_qdqd = derivs.dmod_dvv;
        modID_qdq  = derivs.dmod_dqv';

        dmod_dq_mid2nd  = derivs.dmod_dq;
        dmod_dqd_mid2nd = derivs.dmod_dv;

        dtau_dq_mid2nd  = derivs.dtau_dq;
        dtau_dqd_mid2nd = derivs.dtau_dv;
    
        
        [~, ~, ~, modFD_qq, modFD_qdqd, modFD_qdq, modFD_tauq] = modFD_second_derivatives( model, q, qd, tau, lambda );
        
        modFD_qq_cs   = complexStepJacobian( @(x) outputSelect(1,@modFD_derivatives,model,x,qd,tau,lambda),q );
        modFD_qdqd_cs = complexStepJacobian( @(x) outputSelect(2,@modFD_derivatives,model,q,x,tau,lambda),qd );
        modFD_qdq_cs  = complexStepJacobian( @(x) outputSelect(2,@modFD_derivatives,model,x,qd,tau,lambda),q );
        modFD_tauq_cs = complexStepJacobian( @(x) outputSelect(3,@modFD_derivatives,model,x,qd,tau,lambda),q );

        
        checkValue('ID_q'   , dtau_dq      , dtau_dq_cs            ); % Partials of ID w.r.t. q
        checkValue('ID_qd'  , dtau_dqd     , dtau_dqd_cs           ); % Partials of ID w.r.t. qd
        
        checkValue('FD_q'   , dqdd_dq      , dqdd_dq_cs            ); % Partials of FD w.r.t. q
        checkValue('FD_qd'  , dqdd_dqd     , dqdd_dqd_cs           ); % Partials of FD w.r.t. qd
        checkValue('FD_tau'   , dqdd_dtau    , dqdd_dtau_cs          ); % Partials of FD w.r.t. tau
        
        checkValue('H_q'     , Hpartial     , dH_dq                 ); % Partials of H w.r.t. q
        
        disp('====================================');
        
        checkValue('modID_q'   , dmodID_dq      , dmodID_dq_cs            ); % Partials of modID w.r.t. q
        checkValue('modID_qd'  , dmodID_dqd     , dmodID_dqd_cs           ); % Partials of modID w.r.t. qd
        
        checkValue('modFD_q'   , dmodFD_dq      , dmodFD_dq_cs            ); % Partials of modFD w.r.t. q
        checkValue('modFD_qd'  , dmodFD_dqd     , dmodFD_dqd_cs           ); % Partials of modFD w.r.t. qd
        checkValue('modFD_tau' , dmodFD_dtau    , dmodFD_dtau_cs          ); % Partials of modFD w.r.t. qd
        
        checkValue('modID_qq'    , modID_qq      , modID_qq_cs            ); % SO Partials of modID w.r.t. q,q
        checkValue('modID_qdq'   , modID_qdq     , modID_qdq_cs           ); % SO Partials of modID w.r.t. qd,q
        checkValue('modID_qdqd'  , modID_qdqd    , modID_qdqd_cs          ); % SO Partials of modID w.r.t. qd,qd

        checkValue('modFD_qq'    , modFD_qq      , modFD_qq_cs            ); % SO Partials of modID w.r.t. q,q
        checkValue('modFD_qdq'   , modFD_qdq     , modFD_qdq_cs           ); % SO Partials of modID w.r.t. qd,q
        checkValue('modFD_qdqd'  , modFD_qdqd    , modFD_qdqd_cs          ); % SO Partials of modID w.r.t. qd,qd
        checkValue('modFD_tauq'  , modFD_tauq    , modFD_tauq_cs           ); % SO Partials of modID w.r.t. tau,q
 
        
        disp('====================================');
        checkValue('mod_q'   , dmod_dq_mid2nd    , dmodID_dq_cs           ); % SO Partials of modID w.r.t. qd,qd
        checkValue('mod_qd'  , dmod_dqd_mid2nd   , dmodID_dqd_cs          ); % SO Partials of modID w.r.t. qd,qd

        checkValue('ID_q'   , dtau_dq_mid2nd    , dtau_dq_cs           ); % SO Partials of modID w.r.t. qd,qd
        checkValue('ID_qd'  , dtau_dqd_mid2nd   , dtau_dqd_cs          ); % SO Partials of modID w.r.t. qd,qd
   
    
    end
    
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