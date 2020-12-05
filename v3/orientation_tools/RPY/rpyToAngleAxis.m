function [ angle_axis ] = rpyToAngleAxis( rpy )
angle_axis = rotToAngleAxis(rpyToRot(rpy));
end

