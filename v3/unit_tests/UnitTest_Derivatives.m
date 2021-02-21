% Unit tests for derivatives functions

N = 6;

% Create a random model with N links
model = autoTree(N, 1, pi/3);
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

    % Calculate dynamics quanitites
    [tau]      = ID(model, q ,qd ,qdd);                    % Inverse dynamics
    qdd_ABA    = FDab( model, q, qd, tau);                 % Forward dyanmics
    [H, Cqd]   = HandC(model, q, qd);                      % Mass matrix and bias
    [~, tau_g] = HandC(model, q, qd*0);                    % Gravity force
    Cqd        = Cqd-tau_g;                                % Coriolis force
    [C,Hdot,H2]= CoriolisMatrix( model, q, qd);            % Coriolis matrix
    
    checkValue('H'     , H      , H2                    ); % Mass matrix output is correct
    checkValue('Cqd'   , C*qd   , Cqd                   ); % Generalized Coriolis force
    checkValue('Hdot'  , Hdot   , C+C'                  ); % Hdot -2C skew symmetric
    checkValue('qdd'   , qdd    , qdd_ABA               ); % Forward dynamics
    checkValue('tau'   , tau    , H*qdd+C*qd+tau_g      ); % Inverse dynamics

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
        
        % Check dtau_dq, dtau_dqd, dH_dq
        for i =1:model.NV
            qd_cs = qd;
            q_cs = q;
            qd_cs(i) = qd_cs(i) + sqrt(-1)*eps;
            q_cs(i)  = q(i) + sqrt(-1)*eps;
            
            tau_qd_cs = ID(model, q ,qd_cs ,qdd);
            tau_q_cs  = ID(model, q_cs ,qd ,qdd);
            
            [H_cs, ~] = HandC(model, q_cs, qd);
            
            dH_dq(:,:,i) = imag(H_cs)/eps;
            dtau_dqd_cs(:,i) = imag(tau_qd_cs)/eps;
            dtau_dq_cs(:,i) = imag(tau_q_cs)/eps;
        end
        
        [dtau_dq, dtau_dqd] = ID_derivatives( model, q, qd, qdd );
        
        checkValue('dtau_dq'   , dtau_dq      , dtau_dq_cs            ); % Partials of ID w.r.t. q
        checkValue('dtau_dqd'  , dtau_dqd     , dtau_dqd_cs           ); % Partials of ID w.r.t. qd
        checkValue('dH_dq'     , Hpartial     , dH_dq                 ); % Partials of H w.r.t. q

    end
    
   
    fprintf('\n');
end

function checkValue(name, v1, v2, tolerance)
    if nargin == 3
        tolerance = sqrt(eps);
    end
    value = norm(v1(:)-v2(:));
    fprintf('%10s \t %e\n',name,value);
%     if value > tolerance
%         error('%s is out of tolerance',name);
%     end
end