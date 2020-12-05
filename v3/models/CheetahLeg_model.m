function robot = CheetahLeg_model()
%% robot = CheetahLeg_model()
% Creates a robot model of the Cheetah 3 that is compatible with the
% spatial_v2 dynamics library.
       
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

robot.gravity = [9.81 0 0];
                                        
lo = 0.044;
l1 = .342;

aa3 = sym('a3','real');
dd3 = sym('d3','real');


       
% Modified Denavit Hartenberg parameters alpha_{i-1}, a_{i-1}, and d_i 
% as in J.J.Craig's classic text
dhTable = [0    0   0   ; 
           pi/2 0   lo  ; 
           0    l1  0   ];
       
       
dhTable_sym = [0    0   0   ; 
           pi/2 0   dd3  ; 
           0    aa3  0   ];
       
       
robot.NB = 3;
robot.parent = [0 1 2];
robot.jtype = {'Rz', 'Rz', 'Rz'};
robot.jtype_rotor = {'Rz', 'Rz', 'Rz'};
robot.has_rotor = [1 1 1];


for i = 1:3
    robot.gr{i} = 10;
    alpha = dhTable(i,1);
    a = dhTable(i,2);
    d = dhTable(i,3);
    
    robot.Xtree{i} = xlt([0 0 d]) * xlt([a 0 0]) * rotx(alpha);
    robot.Xtree_motor{i} = xlt([0 0 d]) * xlt([a 0 0]) * rotx(alpha);
    robot.motor_constraint{i} = P;
    
    
    alpha = dhTable(i,1);
    a_s = dhTable_sym(i,2);
    d_s = dhTable_sym(i,3);
    robot.Xtree_sym{i} = xlt([0 0 d_s]) * xlt([a_s 0 0]) * round(rotx(alpha));
    robot.Xtree_motor_sym{i} = robot.Xtree_sym{i};
    
end

% Knee motor sits at the hip
robot.Xtree_motor{3} = eye(6);
robot.I = { eye(6),eye(6),eye(6)};
robot.I_rotor = { eye(6),eye(6),eye(6)};