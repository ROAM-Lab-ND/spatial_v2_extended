% Sanity Check Example for Coriolis Matrix Algorithm

for N = [10 20 30]
    % Create a random model with N links
    model = autoTree(N, 1.5, pi/3);
    model = postProcessModel(model);
    checkDynamics(model,sprintf('N=%d',N));
end

function checkDynamics(model, desc)
    fprintf('====================================\n');
    fprintf('%s\n',desc);
    fprintf('====================================\n');
    model.gravity = [0 0 0]'; 
    
    eCqd  = 0;
    eHdot = 0;
    eGamma = 0;
    
    for i = 1:100
        % Random inertial properties
        for i = 1:model.NB
            model.I{i} = inertiaVecToMat( rand(10,1) );
        end
    
        % Random configuration and velocity
        q   = rand(model.NQ,1)*2*pi;
        q   = normalizeConfVec(model, q); 

        qd  = rand(model.NV,1)*10;
        qdd = zeros(model.NV,1);

        % Calculate dynamics quanitites
        Cqd        = ID(model, q ,qd ,qdd);                    % Inverse dynamics
        [C,Hdot]   = CoriolisMatrix( model, q, qd);            % Coriolis matrix

        eCqd = max(eCqd  , norm(Cqd-C*qd,Inf) );
        eHdot= max(eHdot , norm(Hdot-C-C',Inf) );
        
        % Check Christoffel
        if ~any(model.nv > 1) && ~any(model.has_rotor)
            Gamma = Christoffel(model,q);
            C2 = 0*C;
            for i = 1:model.NB
                C2 = C2 + Gamma(:,:,i)*qd(i);
            end
            eGamma = max(eGamma,norm(C2-C,Inf));
        end
    end
    
    checkValue('Cqd'   , eHdot  , 0                  ); % Generalized Coriolis force
    checkValue('Hdot'  , eGamma , 0                  ); % Hdot -2C skew symmetric
    checkValue('Gamma' , eCqd   , 0                  ); % Christoffel
    
    
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