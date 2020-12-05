function cayley = angleAxisToCayley(ang_axis)

[angle, axis] = angleAxisToAngleAndAxis(ang_axis);

cayley= tan(angle/2)*axis;

end

