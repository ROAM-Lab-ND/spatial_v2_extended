function [model, graphics] = Arm6LinkModel()

% This function creates both the rigid body model struct and a graphics
% cell array that specifies the full teleop arm model

% This version has no rotors and no absolute triplet.

% The 0-configuration for the robot is completely vertical with all joint
% frames / body frames aligned
% The inertial coordinates have +z up (i.e. gravity = [0 0 -9.81 m/s^2])
% The body coordinates have +x forward, +y left, +z up

% remove ground contact body and points initialization

%% Nominal arm parameters

% masses (approx from URDF)
link_masses = [2.93, 1.74, 4.99, 0.82, 0.17, 0.65]; % added 600g for gripper to last mass

% link lengths/offsets (from URDF) (first is at origin, then from prev. joint origin)
joint_offsets = [0, 0, 0; %0, 0, 0.051;
                 0, 0, 0.106;
                 0, 0, 0.071;
                 0, -0.0095, 0.3855;
                 0, 0, 0.362;
                 0, 0.004, 0.035];

% link centers (approx from CAD) used for bounding boxes, inertia estimates
% TODO: improve these, using joint_offsets/2 for now
link_centers = [joint_offsets(2,:)/2;
                joint_offsets(3,:)/2;
                joint_offsets(4,:)/2;
                joint_offsets(5,:)/2;
                joint_offsets(6,:)/2;
                0, 0, 0.100];

% bounding boxes (approx from CAD)
base_box = [0.125, 0.2, 0.175]; % all box lengths are [x,y,z] % TODO: add center locations here?
shoulder_box = [0.130, 0.138, 0.140];
upper_arm_box = [0.130, 0.140, 0.430];
lower_arm_box = [0.070, 0.050, 0.420];
wrist_box = [0.040, 0.060, 0.060];
gripper_box = [0.050, 0.110, 0.200]; % represents wrist roll output

%% Initialize model

% There are 6 links/bodies and 6 rotors
Nb = 0; % start body index at 0
Nl = 0;
model = RBD_model();

%% Shoulder 1 (Shoulder yaw output link)

Nb = Nb+1;
Nl = Nl+1;
model.parent(Nb) = Nb-1;

model.joint{Nb} = revoluteJoint();
model.joint{Nb}.jointAxis = 'z';

model.Xtree{Nb}(1:6,:) = plux(eye(3), joint_offsets(1,:)');
model.I{Nb} = mcI(link_masses(1),link_centers(1,:),boxInertia(link_masses(1),base_box));

model.XtreeKin{Nl} = model.Xtree{Nb};
model.Il{Nl} = model.I{Nb};

model.I_RB{Nb} = model.I{Nb};

graphics{Nb}.boundCenter = link_centers(1,:)';
graphics{Nb}.boundAxes   = base_box/2*1.8;

%% Shoulder 2 (Shoulder pitch output link)

Nb = Nb+1;
Nl = Nl+1;
model.parent(Nb) = Nb-1;

model.joint{Nb} = revoluteJoint();
model.joint{Nb}.jointAxis = 'y';

model.Xtree{Nb}(1:6,:) = plux(eye(3), joint_offsets(2,:)');
model.I{Nb} = mcI(link_masses(2),link_centers(2,:),boxInertia(link_masses(2),shoulder_box));

model.XtreeKin{Nl} = model.Xtree{Nb};
model.Il{Nl} = model.I{Nb};

model.I_RB{Nb} = model.I{Nb};

graphics{Nb}.boundCenter = link_centers(2,:)';
graphics{Nb}.boundAxes   = shoulder_box/2*1.8;

%% Upper arm (Shoulder roll output link)

Nb = Nb+1;
Nl = Nl+1;
model.parent(Nb) = Nb-1;

model.joint{Nb} = revoluteJoint();
model.joint{Nb}.jointAxis = 'z';

model.Xtree{Nb}(1:6,:) = plux(eye(3), joint_offsets(3,:)');
model.I{Nb} = mcI(link_masses(3),link_centers(3,:),boxInertia(link_masses(3),upper_arm_box));

model.XtreeKin{Nl} = model.Xtree{Nb};
model.Il{Nl} = model.I{Nb};

model.I_RB{Nb} = model.I{Nb};

graphics{Nb}.boundCenter = link_centers(3,:)';
graphics{Nb}.boundAxes   = upper_arm_box/2*1.8;

%% Lower arm (Elbow output link)

Nb = Nb+1;
Nl = Nl+1;
model.parent(Nb) = Nb-1;

model.joint{Nb} = revoluteJoint();
model.joint{Nb}.jointAxis = 'y';

model.Xtree{Nb}(1:6,:) = plux(eye(3), joint_offsets(4,:)');
model.I{Nb} = mcI(link_masses(4),link_centers(4,:),boxInertia(link_masses(4),lower_arm_box));

model.XtreeKin{Nl} = model.Xtree{Nb};
model.Il{Nl} = model.I{Nb};

model.I_RB{Nb} = model.I{Nb};

graphics{Nb}.boundCenter = link_centers(4,:)';
graphics{Nb}.boundAxes   = lower_arm_box/2*1.8;

%% Wrist (Wrist pitch output link)

Nb = Nb+1;
Nl = Nl+1;
model.parent(Nb) = Nb-1;

model.joint{Nb} = revoluteJoint();
model.joint{Nb}.jointAxis = 'y';

model.Xtree{Nb}(1:6,:) = plux(eye(3), joint_offsets(5,:)');
model.I{Nb} = mcI(link_masses(5),link_centers(5,:),boxInertia(link_masses(5),wrist_box));

model.XtreeKin{Nl} = model.Xtree{Nb};
model.Il{Nl} = model.I{Nb};

model.I_RB{Nb} = model.I{Nb};

graphics{Nb}.boundCenter = link_centers(5,:)';
graphics{Nb}.boundAxes   = wrist_box/2*1.8;

%% Gripper (Wrist roll output link)

Nb = Nb+1;
Nl = Nl+1;
model.parent(Nb) = Nb-1;

model.joint{Nb} = revoluteJoint();
model.joint{Nb}.jointAxis = 'z';

model.Xtree{Nb}(1:6,:) = plux(eye(3), joint_offsets(6,:)');
model.I{Nb} = mcI(link_masses(6),link_centers(6,:),boxInertia(link_masses(6),gripper_box));

model.XtreeKin{Nl} = model.Xtree{Nb};
model.Il{Nl} = model.I{Nb};

model.I_RB{Nb} = model.I{Nb};

graphics{Nb}.boundCenter = link_centers(6,:)';
graphics{Nb}.boundAxes   = gripper_box/2*1.8;

%% Wrap up
model.NB = Nb;
model.NL = Nl;
model = model.postProcessModel();

end

function I = boxInertia(mass, x) % inertia for modeling as box of mass m and lengths x with uniform density
    I = (norm(x)^2*eye(3) - diag(x.^2))*mass/12;
end