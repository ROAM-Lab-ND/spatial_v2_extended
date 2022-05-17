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
    
    model.joint{i}.gearRatio{1} = 5; % first rotor
    model.joint{i}.gearRatio{2} = 2.4; % second rotor
    model.joint{i}.beltRatio = 1.5;
    
end
model = model.postProcessModel();

q = rand(model.NQ,1); % using q's for the joints as relative link angles

q = model.normalizeConfVec(q);
qd = rand(model.NV,1);
qdd = rand(model.NV,1);

tau = ID(model,q,qd,qdd); % tau gets weird. See notes.

qdd2 = FDab(model,q,qd,tau);
eqdd = norm(qdd-qdd2) % Check consistency of ABA / ID

model.gravity = [0;0;0];
Cqd = ID(model,q,qd,qdd*0);
[C Hdot H] = CoriolisMatrix(model,q,qd);

[H2, Cqd2 ] = HandC(model,q,qd);

eH = norm(H-H2)

eC = norm( Cqd-C*qd ) % Check coriolis matrix
eC2= norm(Cqd2 - Cqd)

eSkew = norm(Hdot-C-C')

Hqdd = ID(model,q,qd*0,qdd); 
eHqdd = norm(H*qdd-Hqdd) % Check mass matrix