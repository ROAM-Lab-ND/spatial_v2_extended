function quat = angleAxisToQuat(ang_axis)

sinc2 = @(z) sinc(z/pi);

%[angle, axis] = angleAxisToAngleAndAxis(ang_axis);
angle = sqrt(ang_axis'*ang_axis);
quat = [cos(angle/2) ; sinc2(angle/2)*ang_axis/2];

end

