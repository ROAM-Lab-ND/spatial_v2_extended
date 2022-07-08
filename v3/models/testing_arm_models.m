clear

%% Test without rotors
disp('Regular old revolutes')
org_model = autoTree(8,2,.3);

Xtree = org_model.Xtree;

model =  RBD_model();
jointAxes = {'x','y','z','x','y','z','x','y','z','x','y','z'};

model.NB = org_model.NB;
model.parent = org_model.parent;


model.Xtree = org_model.Xtree;
model.I     = org_model.I;
    

model.joint{1} = floatingBaseJoint();
for i = 2:model.NB
    model.joint{i} = revoluteJoint();
    model.joint{i}.axis = jointAxes{i}; 
end
model = model.postProcessModel();

q = rand(model.NQ,1);
q = model.normalizeConfVec(q);
qd = rand(model.NV,1);
qdd = rand(model.NV,1);

tau = ID(model,q,qd,qdd);
qdd2 = FDab(model,q,qd,tau);
eqdd = norm( qdd-qdd2 )

model.gravity = [0;0;0];
Cqd = ID(model,q,qd,qdd*0);
[C Hdot H] = CoriolisMatrix(model,q,qd);
eC = norm( Cqd-C*qd )
eSkew = norm(Hdot-C-C')

Hqdd = ID(model,q,qd*0,qdd);
eH = norm(H*qdd-Hqdd)


%% Test without rotors
disp('Absolute Pair Joints')
org_model = autoTree(8,2,.3);

Xtree = org_model.Xtree;

model =  RBD_model();
jointAxes = {'x','y','z','x','y','z','x','y','z','x','y','z'};

model.NB = org_model.NB;
model.parent = org_model.parent;


model.Xtree = org_model.Xtree;
model.I     = org_model.I;
    

model.joint{1} = floatingBaseJoint();
for i = 2:model.NB
    model.joint{i} = revolutePairAbsolute();
    model.joint{i}.axis{1} = jointAxes{i}; 
    model.joint{i}.axis{2} = jointAxes{i}; 
    model.joint{i}.XtreeInternal = randomAdSE3();
    model.I{i} = [model.I{i} zeros(6); zeros(6) model.I{i}];
end
model = model.postProcessModel();

q = rand(model.NQ,1);
q = model.normalizeConfVec(q);
qd = rand(model.NV,1);
qdd = rand(model.NV,1);

tau = ID(model,q,qd,qdd);
qdd2 = FDab(model,q,qd,tau);
eqdd = norm(qdd-qdd2)

model.gravity = [0;0;0];
Cqd = ID(model,q,qd,qdd*0);
[C Hdot H] = CoriolisMatrix(model,q,qd);
eC = norm( Cqd-C*qd )
eSkew = norm(Hdot-C-C')

Hqdd = ID(model,q,qd*0,qdd);
eH = norm(H*qdd-Hqdd)




%% Test with rotors
disp('Rotors')
clear
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
    model.joint{i} = revoluteJointWithRotor();
    model.joint{i}.jointAxis = jointAxes{i};
    model.joint{i}.rotorAxis = jointAxes{i};
    
    model.Xtree{i}(1:6,:) = Xtree{i};
    model.Xtree{i}(7:12,:) = org_model.Xrotor{i};
    model.joint{i}.gearRatio = 5;
    model.I{i} = [model.I{i} zeros(6) ; zeros(6) org_model.I_rotor{i}];
end

model = model.postProcessModel();

q = rand(model.NQ,1);
q = model.normalizeConfVec(q);
qd = rand(model.NV,1);
qdd = rand(model.NV,1);

tau = ID(model,q,qd,qdd);
qdd2 = FDab(model,q,qd,tau);
eqdd = norm( qdd-qdd2 )

model.gravity = [0;0;0];
Cqd = ID(model,q,qd,qdd*0);
[C Hdot H] = CoriolisMatrix(model,q,qd);
eC = norm( Cqd-C*qd )
eSkew = norm(Hdot-C-C')

Hqdd = ID(model,q,qd*0,qdd);
eH = norm(H*qdd-Hqdd)

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
    model.joint{i}.jointAxis{1} = jointAxes{i}; 
    model.joint{i}.jointAxis{2} = jointAxes{i}; 
    
    model.joint{i}.rotorAxis{1} = jointAxes{i}; 
    model.joint{i}.rotorAxis{2} = jointAxes{i};
    
    model.joint{i}.XtreeInternal = randomAdSE3();
    
    model.Xtree{i}(1:6,:) = Xtree{i};
    model.Xtree{i}(7:12,:) = org_model.Xrotor{i};
    model.Xtree{i}(13:18,:) = randomAdSE3();
    
    model.I{i}(1:6,1:6) = randomInertia();
    model.I{i}(7:12,7:12) = randomInertia();
    model.I{i}(13:18,13:18) = randomInertia();
    model.I{i}(19:24,19:24) = randomInertia();
    
    model.joint{i}.gearRatio{1} = 5;
    model.joint{i}.gearRatio{2} = 2.4;
    model.joint{i}.beltRatio{1} = 1.5;
    model.joint{i}.beltRatio{2} = 2.2;
    
end
model = model.postProcessModel();

q = rand(model.NQ,1);
q = model.normalizeConfVec(q);
qd = rand(model.NV,1);
qdd = rand(model.NV,1);

tau = ID(model,q,qd,qdd);
qdd2 = FDab(model,q,qd,tau);
eqdd = norm(qdd-qdd2)

model.gravity = [0;0;0];
Cqd = ID(model,q,qd,qdd*0);
[C Hdot H] = CoriolisMatrix(model,q,qd);
eC = norm( Cqd-C*qd )
eSkew = norm(Hdot-C-C')

Hqdd = ID(model,q,qd*0,qdd);
eH = norm(H*qdd-Hqdd)







