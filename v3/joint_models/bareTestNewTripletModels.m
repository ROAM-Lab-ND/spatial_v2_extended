%% Test with rotors
disp('Absolute Triplet Joints with Rotors')

org_model = autoTree_rotor(8,2,.3); % TODO: up this to 9 bodies?

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
    model.joint{i} = revoluteTripletWithRotors();
    model.joint{i}.jointAxis{1} = jointAxes{i};  % x, y, or z
    model.joint{i}.jointAxis{2} = jointAxes{i};  % x, y, or z
    model.joint{i}.jointAxis{3} = jointAxes{i};  % x, y, or z

    model.joint{i}.rotorAxis{1} = jointAxes{i}; % x, y, or z
    model.joint{i}.rotorAxis{2} = jointAxes{i}; % x, y, or z
    model.joint{i}.rotorAxis{3} = jointAxes{i}; % x, y, or z

    model.joint{i}.XtreeInternal(1:6,:) = randomAdSE3(); % from first link to second link
    model.joint{i}.XtreeInternal(7:12,:) = randomAdSE3(); % from second link to third link

    model.Xtree{i}(1:6,:) = Xtree{i}; % from predecessor to first link
    model.Xtree{i}(7:12,:) = org_model.Xrotor{i}; % from predecessor to first rotor
    model.Xtree{i}(13:18,:) = randomAdSE3(); % from predecessor to second rotor
    model.Xtree{i}(19:24,:) = randomAdSE3(); % from predecessor to third rotor

    model.I{i}(1:6,1:6) = randomInertia(); % first link
    model.I{i}(7:12,7:12) = randomInertia(); % first rotor
    model.I{i}(13:18,13:18) = randomInertia(); % second link
    model.I{i}(19:24,19:24) = randomInertia(); % second rotor
    model.I{i}(25:30,25:30) = randomInertia(); % third rotor
    model.I{i}(31:36,31:36) = randomInertia(); % third link
    
    a(end+1 : end+10) = inertiaMatToVec(model.I{i}(1:6,1:6));
    a(end+1 : end+10) = inertiaMatToVec(model.I{i}(7:12,7:12));
    a(end+1 : end+10) = inertiaMatToVec(model.I{i}(13:18,13:18));
    a(end+1 : end+10) = inertiaMatToVec(model.I{i}(19:24,19:24));
    a(end+1 : end+10) = inertiaMatToVec(model.I{i}(25:30,25:30));
    a(end+1 : end+10) = inertiaMatToVec(model.I{i}(31:36,31:36));

    model.joint{i}.gearRatio{1} = rand()*5; % first rotor
    model.joint{i}.gearRatio{2} = rand()*5; % second rotor
    model.joint{i}.gearRatio{3} = rand()*5; % third rotor
    model.joint{i}.beltRatio{1} = rand()*5;
    model.joint{i}.beltRatio{2} = rand()*5;
    model.joint{i}.beltRatio{3} = rand()*5;

end

model = model.postProcessModel();

q = rand(model.NQ,1); % using q's for the joints as relative link angles

q = model.normalizeConfVec(q);
qd = rand(model.NV,1);
qdr = rand(model.NV,1);
qdd = rand(model.NV,1);

tau = ID(model,q,qd,qdd); % tau gets weird. See notes.

qdd2 = FDab(model,q,qd,tau);
eqdd = norm(qdd-qdd2) % Check consistency of ABA / ID

Y = RegressorClassical(model, q,qd,qdd);
eY = norm( tau - Y*a)


Y_sl = RegressorSL(model, q, qd, qdr, qdd);
tau_sl = ID_SlotineLi(model, q, qd, qdr,qdd);

eY_SL = norm( tau_sl - Y_sl*a)


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

Hinv = Hinverse(model,q);
eHinv = norm(Hinv - inv(H))







