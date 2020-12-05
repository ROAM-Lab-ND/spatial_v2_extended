function mrp = angleAxisToMRP(ang_axis)

[angle, axis] = angleAxisToAngleAndAxis(ang_axis);
mrp= tan(angle/4)*axis;

end

