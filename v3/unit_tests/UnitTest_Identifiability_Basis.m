%  This is intended to verify the matrices (Minimal_Basis and Perp_Basis)  obtained from 
% "RPNA_URDF_Examples.m".
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

% panda basis from DH
[N1, M1, V1, C1] = RPNA(model_DH, 0);
[~, Minimal_Basis_DH, Perp_Basis_DH] = ComputeBases(model_DH, N1, M1); % After this step, the perp basis is correct
Perp_Basis_DH = rref(Perp_Basis_DH')'; % However, if you want to use it to project to minimal parameters, you need this step

% panda basis from URDF
[N2, M2, V2, C2] = RPNA(model_URDF, 0);
[~, Minimal_Basis_URDF, Perp_Basis_URDF] = ComputeBases(model_URDF, N2, M2); % After this step, the perp basis is correct
Perp_Basis_URDF = rref(Perp_Basis_URDF')'; % However, if you want to use it to project to minimal parameters, you need this step

%% Verification loop
timeDuration = 2.5;
currentTime = 0;
timeStep = 0.005;
for i =0:timeDuration/timeStep 
    q = -cos(currentTime)*ones(dof,1);
    qd = sin(currentTime)*ones(dof,1);
    qdd = cos(currentTime)*ones(dof,1);
    qd_r = qd;
    disp("===============================================")
    
    tau1 = inverseDynamics(robotics_toolbox_robot, q, qd, qdd);
    
    [tau2, ~] = ID( model_DH, q, qd, qdd );
    checkValue('ID with model from DH',tau2, tau1, 1e-5);

    [tau3, ~] = ID( model_URDF, q, qd, qdd );
    checkValue('ID with model from URDF',tau3, tau1, 1e-5);

    [Y_DH, ~] = RegressorClassical( model_DH, q, qd,qdd);
    tau4 = Y_DH*W;
    checkValue('Torque of model from DH',tau4, tau1, 1e-5);
    
    tau5 = Y_DH * Minimal_Basis_DH * Perp_Basis_DH' *W;
    checkValue('Torque of model from DH with basis',tau5, tau1, 1e-5);

    [Y_URDF, ~] = RegressorClassical( model_URDF, q, qd,qdd);
    tau6 = Y_URDF*W;
    checkValue('Torque of model from URDF',tau6, tau1, 1e-5);
    
    tau7 = Y_URDF * Minimal_Basis_URDF * Perp_Basis_URDF' *W;
    checkValue('Torque of model from URDF with basis',tau7, tau1, 1e-5);

    currentTime = currentTime + timeStep;
end

function checkValue(name, v1, v2, tolerance)
    if nargin == 3
        tolerance = sqrt(eps);
    end
    value = norm(v1(:)-v2(:));
    fprintf('%s \t %e\n',name,value);
    if value > tolerance
        error('%s is out of tolerance',name);
    end
end