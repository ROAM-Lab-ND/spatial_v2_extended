function [model, graphics] = Arm6LinkAbsModel()

% This function creates both the rigid body model struct and a graphics
% cell array that specifies the full teleop arm model

% This version has rotors and the absolute triplet for the final three
% joints.

% The 0-configuration for the robot is completely vertical with all joint
% frames / body frames aligned
% The inertial coordinates have +z up (i.e. gravity = [0 0 -9.81 m/s^2])
% The body coordinates have +x forward, +y left, +z up

% remove ground contact body and points initialization

%% Nominal arm parameters

% motor sizes
U12_dia = 0.120;
U12_width = 0.04;

U10_dia = 0.100;
U10_width = 0.035;

% masses (approx from URDF)
link_masses = [2.93, 1.74, 4.99, 0.82, 0.17, 0.65]; % added 600g for gripper to last mass
rotor_masses = [0.1, 0.06]; % TODO: update these from CAD? U10 rotor is small on purpose for now

% link lengths/offsets (from URDF) (first is at origin, then from prev. joint origin)
joint_offsets = [0, 0, 0; %0, 0, 0.051;
                 0, 0, 0.106;
                 0, 0, 0.071;
                 0, -0.0095, 0.3855;
                 0, 0, 0.362;
                 0, 0.004, 0.035];

% rotor offsets (from CAD) (from corresponding joint origin for joints 1-4, from upper arm origin for 5 and 6)
rotor_offsets = [0, 0, -0.02125;
                0, 0.09625, 0;
                0, 0, -0.02125;
                0, 0.07088, -0.307;
                0, -0.01, -0.307; % relative to elbow frame
                0, 0.01, -0.205]; % relative to elbow frame

rotor_offsets(1:4,:) = rotor_offsets(1:4,:) + joint_offsets(1:4,:);
rotor_offsets(5,:) = rotor_offsets(5,:) + joint_offsets(4,:);
rotor_offsets(6,:) = rotor_offsets(6,:) + joint_offsets(4,:);

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
Nl = 0; % start link index at 0
Nrb = 0;
model = RBD_model();


%% Shoulder 1 (Shoulder yaw output link)

Nb = Nb+1;
Nl = Nl+1;
model.parent(Nb) = Nb-1;

model.joint{Nb} = revoluteJointWithRotor();
model.joint{Nb}.jointAxis = 'z';
model.joint{Nb}.rotorAxis = 'z';

model.Xtree{Nb}(1:6,:) = plux(eye(3), joint_offsets(1,:)');
model.Xtree{Nb}(7:12,:) = plux(eye(3), rotor_offsets(1,:)');
model.joint{Nb}.gearRatio = 6;
model.I{Nb} = [mcI(link_masses(1),link_centers(1,:),boxInertia(link_masses(1),base_box)),  zeros(6); zeros(6), mcI(rotor_masses(1),[0 0 0],boxInertia(rotor_masses(1),[U12_dia, U12_dia, U12_width]))];

model.XtreeKin{Nl} = plux(eye(3), joint_offsets(1,:)');
model.XtreeKinRot{Nl} = plux(eye(3), rotor_offsets(1,:)');
model.Il{Nl} = mcI(link_masses(1),link_centers(1,:),boxInertia(link_masses(1),base_box));
model.Ir{Nl} = mcI(rotor_masses(1),[0 0 0],boxInertia(rotor_masses(1),[U12_dia, U12_dia, U12_width]));

Nrb = Nrb+1;
model.I_RB{Nrb} = model.Il{Nl};
Nrb = Nrb+1;
model.I_RB{Nrb} = model.Ir{Nl};

graphics{Nl}.boundCenter = link_centers(1,:)';
graphics{Nl}.boundAxes   = base_box/2*1.8;
graphics{Nl}.boundCenterRot = [0 0 0]';
graphics{Nl}.boundAxesRot   = ([U12_dia, U12_dia, U12_width]/2)*1.8;

%% Shoulder 2 (Shoulder pitch output link)

Nb = Nb+1;
Nl = Nl+1;
model.parent(Nb) = Nb-1;

model.joint{Nb} = revoluteJointWithRotor();
model.joint{Nb}.jointAxis = 'y';
model.joint{Nb}.rotorAxis = 'y';

model.Xtree{Nb}(1:6,:) = plux(eye(3), joint_offsets(2,:)');
model.Xtree{Nb}(7:12,:) = plux(eye(3), rotor_offsets(2,:)');
model.joint{Nb}.gearRatio = 6;
model.I{Nb} = [mcI(link_masses(2),link_centers(2,:),boxInertia(link_masses(2),shoulder_box)),  zeros(6); zeros(6), mcI(rotor_masses(1),[0 0 0],boxInertia(rotor_masses(1),[U12_dia, U12_width, U12_dia]))];

model.XtreeKin{Nl} = plux(eye(3), joint_offsets(2,:)');
model.XtreeKinRot{Nl} = plux(eye(3), rotor_offsets(2,:)');
model.Il{Nl} = mcI(link_masses(2),link_centers(2,:),boxInertia(link_masses(2),shoulder_box));
model.Ir{Nl} = mcI(rotor_masses(1),[0 0 0],boxInertia(rotor_masses(1),[U12_dia, U12_width, U12_dia]));

Nrb = Nrb+1;
model.I_RB{Nrb} = model.Il{Nl};
Nrb = Nrb+1;
model.I_RB{Nrb} = model.Ir{Nl};

graphics{Nl}.boundCenter = link_centers(2,:)';
graphics{Nl}.boundAxes   = shoulder_box/2*1.8;
graphics{Nl}.boundCenterRot = [0 0 0]';
graphics{Nl}.boundAxesRot   = ([U12_dia, U12_width, U12_dia]/2)*1.8;

%% Upper arm (Shoulder roll output link)

Nb = Nb+1;
Nl = Nl+1;
model.parent(Nb) = Nb-1;

model.joint{Nb} = revoluteJointWithRotor();
model.joint{Nb}.jointAxis = 'z';
model.joint{Nb}.rotorAxis = 'z';

model.Xtree{Nb}(1:6,:) = plux(eye(3), joint_offsets(3,:)');
model.Xtree{Nb}(7:12,:) = plux(eye(3), rotor_offsets(3,:)');
model.joint{Nb}.gearRatio = 6;
model.I{Nb} = [mcI(link_masses(3),link_centers(3,:),boxInertia(link_masses(3),upper_arm_box)),  zeros(6); zeros(6), mcI(rotor_masses(1),[0 0 0],boxInertia(rotor_masses(1),[U12_dia, U12_dia, U12_width]))];

model.XtreeKin{Nl} = plux(eye(3), joint_offsets(3,:)');
model.XtreeKinRot{Nl} = plux(eye(3), rotor_offsets(3,:)');
model.Il{Nl} = mcI(link_masses(3),link_centers(3,:),boxInertia(link_masses(3),upper_arm_box));
model.Ir{Nl} = mcI(rotor_masses(1),[0 0 0],boxInertia(rotor_masses(1),[U12_dia, U12_dia, U12_width]));

Nrb = Nrb+1;
model.I_RB{Nrb} = model.Il{Nl};
Nrb = Nrb+1;
model.I_RB{Nrb} = model.Ir{Nl};

graphics{Nl}.boundCenter = link_centers(3,:)';
graphics{Nl}.boundAxes   = upper_arm_box/2*1.8;
graphics{Nl}.boundCenterRot = [0 0 0]';
graphics{Nl}.boundAxesRot   = ([U12_dia, U12_dia, U12_width]/2)*1.8;

%% Lower arm (Elbow output link)

Nb = Nb+1;
Nl = Nl+1;
model.parent(Nb) = Nb-1;

% model.joint{Nb} = revoluteJointWithRotor();
% model.joint{Nb}.jointAxis = 'y';
% model.joint{Nb}.rotorAxis = 'y';
% 
% model.Xtree{Nb}(1:6,:) = plux(eye(3), joint_offsets(4,:)');
% model.Xtree{Nb}(7:12,:) = plux(eye(3), rotor_offsets(4,:)');
% model.joint{Nb}.gearRatio = 6;
% model.I{Nb} = [mcI(link_masses(4),link_centers(4,:),boxInertia(link_masses(4),lower_arm_box)),  zeros(6); zeros(6), mcI(rotor_masses(1),[0 0 0],boxInertia(rotor_masses(1),[U12_dia, U12_width, U12_dia]))];

% include elbow and wrist pitch here
model.joint{Nb} = revoluteTripletWithRotors();
model.joint{Nb}.jointAxis{1} = 'y';
model.joint{Nb}.rotorAxis{1} = 'y'; 

model.joint{Nb}.jointAxis{2} = 'y'; 
model.joint{Nb}.rotorAxis{2} = 'y';

model.joint{Nb}.jointAxis{3} = 'z'; 
model.joint{Nb}.rotorAxis{3} = 'y';

model.joint{Nb}.XtreeInternal(1:6,:) = plux(eye(3), joint_offsets(5,:)'); % second joint offset
model.joint{Nb}.XtreeInternal(7:12,:) = plux(eye(3), joint_offsets(6,:)'); % third joint offset

model.Xtree{Nb}(1:6,:) = plux(eye(3), joint_offsets(4,:)'); % first joint offset
model.Xtree{Nb}(7:12,:) = plux(eye(3), rotor_offsets(4,:)'); % first rotor offset
model.Xtree{Nb}(13:18,:) = plux(eye(3), rotor_offsets(5,:)'); % second rotor offset
model.Xtree{Nb}(19:24,:) = plux(eye(3), rotor_offsets(6,:)'); % third rotor offset

% inertias: link 1, rotor 1, link 2, rotor 2, rotor 3, link 3
model.I{Nb}(1:6,1:6) = mcI(link_masses(4),link_centers(4,:),boxInertia(link_masses(4),lower_arm_box)); % first link inertia
model.I{Nb}(7:12,7:12) = mcI(rotor_masses(1),[0 0 0],boxInertia(rotor_masses(1),[U12_dia, U12_width, U12_dia])); % first rotor inertia
model.I{Nb}(13:18,13:18) = mcI(link_masses(5),link_centers(5,:),boxInertia(link_masses(5),wrist_box)); % second link inertia
model.I{Nb}(19:24,19:24) = mcI(rotor_masses(2),[0 0 0],boxInertia(rotor_masses(2),[U10_dia, U10_width, U10_dia])); % second rotor inertia 
model.I{Nb}(25:30,25:30) = mcI(rotor_masses(2),[0 0 0],boxInertia(rotor_masses(2),[U10_dia, U10_width, U10_dia])); % third rotor inertia
model.I{Nb}(31:36,31:36) = mcI(link_masses(6),link_centers(6,:),boxInertia(link_masses(6),gripper_box)); % third link inertia

model.joint{Nb}.gearRatio{1} = 6;
model.joint{Nb}.gearRatio{2} = 6;
model.joint{Nb}.gearRatio{3} = 6;
model.joint{Nb}.beltRatio{1} = 1;
model.joint{Nb}.beltRatio{2} = 1;
model.joint{Nb}.beltRatio{3} = 2;

% for elbow output / lower arm link
model.XtreeKin{Nl} = plux(eye(3), joint_offsets(4,:)');
model.XtreeKinRot{Nl} = plux(eye(3), rotor_offsets(4,:)');
model.Il{Nl} = mcI(link_masses(4),link_centers(4,:),boxInertia(link_masses(4),lower_arm_box));
model.Ir{Nl} = mcI(rotor_masses(1),[0 0 0],boxInertia(rotor_masses(1),[U12_dia, U12_width, U12_dia]));

Nrb = Nrb+1;
model.I_RB{Nrb} = model.Il{Nl}; % first link inertia
Nrb = Nrb+1;
model.I_RB{Nrb} = model.Ir{Nl}; % first rotor inertia

graphics{Nl}.boundCenter = link_centers(4,:)';
graphics{Nl}.boundAxes   = lower_arm_box/2*1.8;
graphics{Nl}.boundCenterRot = [0 0 0]';
graphics{Nl}.boundAxesRot   = ([U12_dia, U12_width, U12_dia]/2)*1.8;

Nl = Nl+1; % for wrist pitch output / wrist link
model.XtreeKin{Nl} = plux(eye(3), joint_offsets(5,:)');
model.XtreeKinRot{Nl} = plux(eye(3), rotor_offsets(5,:)'); % Note: this is relative to elbow parent body, not wrist pitch parent body
model.Il{Nl} = mcI(link_masses(5),link_centers(5,:),boxInertia(link_masses(5),wrist_box)); 
model.Ir{Nl} = mcI(rotor_masses(2),[0 0 0],boxInertia(rotor_masses(2),[U10_dia, U10_width, U10_dia]));

Nrb = Nrb+1;
model.I_RB{Nrb} = model.Il{Nl}; % second link inertia
Nrb = Nrb+1;
model.I_RB{Nrb} = model.Ir{Nl}; % second rotor inertia

graphics{Nl}.boundCenter = link_centers(5,:)';
graphics{Nl}.boundAxes   = wrist_box/2*1.8;
graphics{Nl}.boundCenterRot = [0 0 0]';
graphics{Nl}.boundAxesRot  = ([U10_dia, U10_width, U10_dia]/2)*1.8;


Nl = Nl + 1; % for wrist roll output / gripper link
model.XtreeKin{Nl} = plux(eye(3), joint_offsets(6,:)');
model.XtreeKinRot{Nl} = plux(eye(3), rotor_offsets(6,:)');
model.Il{Nl} = mcI(link_masses(6),link_centers(6,:),boxInertia(link_masses(6),gripper_box));
model.Ir{Nl} = mcI(rotor_masses(2),[0 0 0],boxInertia(rotor_masses(2),[U10_dia, U10_width, U10_dia]));

Nrb = Nrb+1;
model.I_RB{Nrb} = model.Ir{Nl}; % third rotor inertia
Nrb = Nrb+1;
model.I_RB{Nrb} = model.Il{Nl}; % third link inertia

graphics{Nl}.boundCenter = link_centers(6,:)';
graphics{Nl}.boundAxes   = gripper_box/2*1.8;
graphics{Nl}.boundCenterRot = [0 0 0]';
graphics{Nl}.boundAxesRot   = ([U10_dia, U10_width, U10_dia]/2)*1.8;

%% Wrist (Wrist pitch output link)

% Nb = Nb+1;
% model.parent(Nb) = Nb-1;

% model.joint{Nb} = revoluteJointWithRotor();
% model.joint{Nb}.jointAxis = 'y';
% model.joint{Nb}.rotorAxis = 'y';
% 
% model.Xtree{Nb}(1:6,:) = plux(eye(3), joint_offsets(5,:)');
% model.Xtree{Nb}(7:12,:) = plux(eye(3), rotor_offsets(5,:)');
% model.joint{Nb}.gearRatio = 6;
% model.I{Nb} = [mcI(link_masses(5),link_centers(5,:),boxInertia(link_masses(5),wrist_box)),  zeros(6); zeros(6), mcI(rotor_masses(2),[0 0 0],boxInertia(rotor_masses(2),[U10_dia, U10_dia, U10_width]))];
% 
% 
% graphics{Nb}.boundCenter = link_centers(5,:)';
% graphics{Nb}.boundAxes   = wrist_box/2*1.8;
% graphics{Nb}.boundCenterRot = [0 0 0]';
% graphics{Nb}.boundAxesRot   = ([U10_dia, U10_dia, U10_width]/2)*1.8;

%% Gripper (Wrist roll output link)

% adding rotor to this one will be extra tricky...need to rotate output
% frame to get positive axis?

% Nb = Nb+1;
% Nl = Nl+1;
% model.parent(Nb) = Nb-1;
% 
% model.joint{Nb} = revoluteJointWithRotor();
% model.joint{Nb}.jointAxis = 'z';
% model.joint{Nb}.rotorAxis = 'y';
% 
% model.Xtree{Nb}(1:6,:) = plux(eye(3), joint_offsets(6,:)');
% model.Xtree{Nb}(7:12,:) = plux(eye(3), rotor_offsets(6,:)');
% model.joint{Nb}.gearRatio = 6*2; % GR is 6, belt is 2
% model.I{Nb} = [mcI(link_masses(6),link_centers(6,:),boxInertia(link_masses(6),gripper_box)),  zeros(6); zeros(6), mcI(rotor_masses(2),[0 0 0],boxInertia(rotor_masses(2),[U10_dia, U10_dia, U10_width]))];
% 
% model.XtreeKin{Nl} = plux(eye(3), joint_offsets(6,:)');
% model.XtreeKinRot{Nl} = plux(eye(3), rotor_offsets(6,:)');
% model.Il{Nl} = mcI(link_masses(6),link_centers(6,:),boxInertia(link_masses(6),gripper_box));
% model.Ir{Nl} = mcI(rotor_masses(2),[0 0 0],boxInertia(rotor_masses(2),[U10_dia, U10_dia, U10_width]));
% 
% graphics{Nl}.boundCenter = link_centers(6,:)';
% graphics{Nl}.boundAxes   = gripper_box/2*1.8;
% graphics{Nl}.boundCenterRot = [0 0 0]';
% graphics{Nl}.boundAxesRot   = ([U10_dia, U10_dia, U10_width]/2)*1.8;


%% Wrap up

model.NB = Nb;
model.NL = Nl;
model = model.postProcessModel();

end

function I = boxInertia(mass, x) % inertia for modeling as box of mass m and lengths x with uniform density
    I = (norm(x)^2*eye(3) - diag(x.^2))*mass/12;
end