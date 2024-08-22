% Joint	    a(link length)     d(link offset)     alpha(link twist)    theta(joint angle)
% Joint 1	    0	                      0.333	              0	                       theta_1
% Joint 2	    0	                      0	                  -pi/2	                   theta_2
% Joint 3	    0	                      0.316	              pi/2	                   theta_3
% Joint 4	    0.0825	              0	                  pi/2	                   theta_4
% Joint 5	    -0.0825	          0.384	              -pi/2	                   theta_5
% Joint 6	    0	                      0	                  pi/2	                   theta_6
% Joint 7	    0.088	              0	                  pi/2	                   theta_7
% Flange	    0	                      0.107	              0	                       0
% https://www.tq-franka.cn/FCI/control_parameters.html

function robot = Panda_model(params)
% Creates a robot model of the Panda that is compatible with the
% spatial_v2 dynamics library.

robot.NB = 7;
robot.gravity = [0 0 -9.81];

% get spatial inertia of each link from the vector of inertia parameters
% a = [m hx hy hz Ixx Iyy Izz Ixy Ixz Iyz]'
if nargin == 1
    for i = 1:robot.NB
        robot.I{i} = inertiaVecToMat(params(10*i-9:10*i));
    end
elseif nargin > 1 
     error("Incorrect number of input arguments!");
end

a4 = 0.0825;
a5 = -0.0825;
a7 = 0.088;
d1 = 0.333;
d3 = 0.316;
d5 = 0.384;
% 
aa4 = sym('a4','real');
aa5 = sym('a5','real');
aa7 = sym('a7','real');
dd1 = sym('d1','real');
dd3 = sym('d3','real');
dd5 = sym('d5','real');

% Modified Denavit Hartenberg parameters alpha_{i-1}, a_{i-1}, and d_i 
% as in J.J.Craig's classic text
dhTable = [0        0   d1   ; 
           -pi/2    0   0   ; 
           pi/2     0   d3  ; 
           pi/2     a4  0  ; 
           -pi/2    a5  d5  ; 
           pi/2     0   0   ;
           pi/2     a7   0 ];
       
       
dhTable_sym = [0        0    dd1   ; 
               -pi/2    0    0   ; 
               pi/2     0    dd3  ; 
               pi/2     aa4  0  ; 
               -pi/2    aa5  dd5  ; 
               pi/2     0    0   ;
               pi/2     aa7  0 ];      

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



robot.parent = [0 1 2 3 4 5 6];
robot.jtype = {'Rz', 'Rz', 'Rz','Rz','Rz','Rz', 'Rz'};

robot.jtype_rotor = {'Rz', 'Rz', 'Rz','Rz','Rz','Rz','Rz'};
robot.has_rotor = zeros(7,1);

for i = 1:7
    robot.gr{i} = 5;
end

robot.Xtree = cell(1,7);
for i = 1:robot.NB
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


