[spatial_v2_model, robotics_toolbox_robot] = URDF_to_spatialv2_model('puma560_robot.urdf', 1);
q = rand(6,1);
qd = rand(6,1);
qdd = rand(6,1);

robotics_toolbox_robot.DataFormat = 'column';
robotics_toolbox_robot.Gravity = [0 0 -9.81]';


M1 = massMatrix(robotics_toolbox_robot, q);
[H, ~] = HandC(spatial_v2_model, q,0*q);

tau1 = inverseDynamics(robotics_toolbox_robot, q, qd, qdd);
tau2 = ID(spatial_v2_model, q, qd, qdd);

checkValue('MassMatrix',M1, H, 1e-5);
checkValue('InverseDyn',tau1,tau2, 1e-5);


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
