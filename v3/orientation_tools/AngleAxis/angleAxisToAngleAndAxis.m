function [angle,axis] = angleAxisToAngleAndAxis(angle_axis)

angle = norm(angle_axis);
if angle == 0
    axis = [0 0 0]';
else
    axis = angle_axis/angle;
end

end

