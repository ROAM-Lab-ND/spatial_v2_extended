function robot = Puma560_model()
%% robot = Puma560_model()
% Creates a robot model of the Puma560 that is compatible with the
% spatial_v2 dynamics library.

a2 = .2;
d3 = .4;
a3 = .153;
d4 = .839;

aa2 = sym('a2','real');
aa3 = sym('a3','real');
dd3 = sym('d3','real');
dd4 = sym('d4','real');


% Modified Denavit Hartenberg parameters alpha_{i-1}, a_{i-1}, and d_i 
% as in J.J.Craig's classic text
dhTable = [0        0   0   ; 
           -pi/2    0   0   ; 
           0        a2  d3  ; 
           -pi/2    a3  d4  ; 
           pi/2     0   0   ; 
           -pi/2    0   0   ];
       
       
dhTable_sym = [0        0   0   ; 
               -pi/2    0   0   ; 
               0        aa2  dd3  ; 
               -pi/2    aa3  dd4  ; 
                pi/2     0   0   ; 
               -pi/2    0   0   ];       

robot.NB = 6;
robot.parent = [0 1 2 3 4 5];
robot.jtype = {'Rz', 'Rz', 'Rz','Rz','Rz','Rz'};

robot.jtype_rotor = {'Rz', 'Rz', 'Rz','Rz','Rz','Rz'};
robot.has_rotor = ones(6,1);

for i = 1:6
    robot.gr{i} = 5;
end


% Constraints on motor parameters to ensure rotational symmetry (about z)
rotor_symmetry_constraint = zeros(6,10);
rotor_symmetry_constraint(1,2) = 1;
rotor_symmetry_constraint(2,3) = 1;
rotor_symmetry_constraint(3,5) = 1; rotor_symmetry_constraint(3,6) = -1;
rotor_symmetry_constraint(4,8) = 1;
rotor_symmetry_constraint(5,9) = 1;
rotor_symmetry_constraint(6,10)= 1; 

B = RangeBasis( null(rotor_symmetry_constraint) ); 
P = B*B'; % Projector onto nullspace of symmetry constraint


robot.Xtree = {  };
for i = 1:6
    alpha = dhTable(i,1);
    a = dhTable(i,2);
    d = dhTable(i,3);
    robot.Xtree{i} = xlt([0 0 d]) * xlt([a 0 0]) * rotx(alpha);
    robot.Xrotor{i} = robot.Xtree{i};
    
    
    alpha = dhTable(i,1);
    a_s = dhTable_sym(i,2);
    d_s = dhTable_sym(i,3);
    robot.Xtree_sym{i} = xlt([0 0 d_s]) * xlt([a_s 0 0]) * round(rotx(alpha));
    robot.Xrotor_sym{i} = robot.Xtree_sym{i};
    
    robot.rotor_constraint{i} = P;
end

robot.gravity = [0 0 -9.81];
robot.I = { eye(6),eye(6),eye(6),eye(6),eye(6),eye(6)};
robot.I_rotor = { eye(6),eye(6),eye(6),eye(6),eye(6),eye(6)};