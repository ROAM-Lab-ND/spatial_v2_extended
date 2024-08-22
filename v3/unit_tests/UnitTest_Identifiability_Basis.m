clear;
clc;
%%
[model_URDF,  robotics_toolbox_robot] = URDF_to_spatialv2_model("panda_arm_no_fixed.urdf");

dof = robotics_toolbox_robot.NumBodies;
W = cell(dof,1);
for i = 1:dof
    Ixx = robotics_toolbox_robot.Bodies{1,i}.Inertia(1);
    Iyy = robotics_toolbox_robot.Bodies{1,i}.Inertia(2);
    Izz = robotics_toolbox_robot.Bodies{1,i}.Inertia(3);
    Iyz = robotics_toolbox_robot.Bodies{1,i}.Inertia(4);
    Ixz = robotics_toolbox_robot.Bodies{1,i}.Inertia(5);
    Ixy = robotics_toolbox_robot.Bodies{1,i}.Inertia(6);
    m = robotics_toolbox_robot.Bodies{1,i}.Mass;
    Cx = robotics_toolbox_robot.Bodies{1,i}.CenterOfMass(1);
    Cy = robotics_toolbox_robot.Bodies{1,i}.CenterOfMass(2);
    Cz = robotics_toolbox_robot.Bodies{1,i}.CenterOfMass(3);
    W{i,1} = [ m , m*Cx, m*Cy, m*Cz, Ixx, Iyy, Izz, Iyz, Ixz, Ixy]';
end
W = cell2mat(W);

model_DH = Panda_model(W);
robotics_toolbox_robot.DataFormat = 'column';
robotics_toolbox_robot.Gravity = [0 0 -9.81]';
%% load basis results of Panda robot modeled by DH or URDF

% from RPNA_Examples
load("panda_basis_from_DH.mat");
Minimal_Basis_DH = Minimal_Basis;
Perp_Basis_DH = Perp_Basis;

% from RPNA_URDF_Examples
load("panda_basis_from_URDF.mat");
Minimal_Basis_URDF = Minimal_Basis;
Perp_Basis_URDF = Perp_Basis;

%% Verification loop
timeDuration = 2.5;
currentTime = 0;
timeStep = 0.005;
for i =0:timeDuration/timeStep 
    q = -cos(currentTime)*ones(dof,1);
    qd = sin(currentTime)*ones(dof,1);
    qdd = cos(currentTime)*ones(dof,1);
    qd_r = qd;
    
    tau1 = inverseDynamics(robotics_toolbox_robot, q, qd, qdd)
    
    [tau2, ~] = ID( model_DH, q, qd, qdd )

    [tau3, ~] = ID( model_URDF, q, qd, qdd )

    [Y_DH, ~] = RegressorClassical( model_DH, q, qd,qdd);
    tau4 = Y_DH*W
    
    tau5 = Y_DH * Minimal_Basis_DH * Perp_Basis_DH' *W

    [Y_URDF, ~] = RegressorClassical( model_URDF, q, qd,qdd);
    tau6 = Y_URDF*W
    
    tau7 = Y_URDF * Minimal_Basis_URDF * Perp_Basis_URDF' *W

    currentTime = currentTime + timeStep;
end