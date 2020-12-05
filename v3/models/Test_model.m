function robot = Test_model()

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



robot.NB = 2;
robot.parent = [0 1];
robot.jtype = {'Rz', 'Rz'};

robot.Xtree = { eye(6), plux(expm(skew([pi/2 0 0])) , [1 0 0]) };

robot.gravity = [0 0 0];		% zero gravity is not the default,
                                        % so it must be stated explicitly
                                        
robot.jtype_rotor = {'Rz', 'Rz'};
robot.has_rotor = ones(2,1);

for i = 1:2
    robot.gr{i} = 5;
    robot.Xrotor{i} = robot.Xtree{i};
    robot.rotor_constraint{i} = P;
    robot.I{i} = eye(6);
    robot.I_rotor{i} = eye(6);
end