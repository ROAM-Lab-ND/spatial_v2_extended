function robot = TX40_model()
%% robot = TX40_model()
% Creates a robot model of the TX40 that is compatible with the
% spatial_v2 dynamics library.

d3 = .225;
r3 = 0.035;
r4 = .225;

% Modified Denavit Hartenberg parameters alpha_{i-1}, a_{i-1}, and d_i 
% as in J.J.Craig's classic text
dhTable = [0        0       0   ; 
           -pi/2    0       0   ; 
           0        d3      r3  ; 
           pi/2     0       r4  ;
           -pi/2    0       0   ;
           pi/2     0       0   ];

% Constraints on motor parameters to ensure rotational symmetry (about z)
motor_symmetry_constraint = zeros(6,10);
motor_symmetry_constraint(1,2) = 1;
motor_symmetry_constraint(2,3) = 1;
motor_symmetry_constraint(3,5) = 1; motor_symmetry_constraint(3,6) = -1;
motor_symmetry_constraint(4,8) = 1;
motor_symmetry_constraint(5,9) = 1;
motor_symmetry_constraint(6,10)= 1; 


B = RangeBasis( null(motor_symmetry_constraint) ); 
P = B*B'; % Projector onto nullspace of symmetry constraint

robot.NB = 6;
robot.parent = [0 1 2 3 4 5];
robot.jtype = {'Rz', 'Rz', 'Rz','Rz','Rz','Rz'};
robot.jtype_rotor = {'Rz', 'Rz', 'Rz','Rz','Rz','Rz'};
robot.has_rotor   = ones(6,1);

robot.gr = { 32,       32,       45 ,     48    ,    45     ,  32  };

robot.Xtree = {  };
for i = 1:6
    alpha = dhTable(i,1);
    a = dhTable(i,2);
    d = dhTable(i,3);
    robot.Xtree{i} = xlt([0 0 d]) * xlt([a 0 0]) * rotx(alpha);
    robot.Xrotor{i} = robot.Xtree{i};
    robot.motor_constraint{i} = P;
end

robot.gravity = [0 0 0];		% zero gravity is not the default,
for i = 1:robot.NB
    robot.I{i} = eye(6);
    robot.I_rot{i} = eye(6);
end
