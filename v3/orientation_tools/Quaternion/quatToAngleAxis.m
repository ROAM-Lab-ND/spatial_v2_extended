function [ angle_axis ] = quatToAngleAxis( quat )

angle_axis = rotToAngleAxis(quatToRot(quat));
 
end

