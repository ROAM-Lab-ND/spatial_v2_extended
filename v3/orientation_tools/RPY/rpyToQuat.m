function [ quat ] = rpyToQuat( rpy )
% rpyToQuat takes earth-fixed sequenial roll-pitch-yaw euler angles to a
% unit quaternion
%   [quat] = rpyToQuat(rpy)
%   quat= [q0,q1,q2,q3] assumed to follow q0 = cos(angle / 2) while
%   [q1,q2,q3] = axis*sin(angle / 2) for angle/axis expression of R
%   note: earth-fixed roll-pitch-yaw is the same as body-fixed
%   yaw-pitch-roll sequence 

    roll = rpy(1);
    pitch = rpy(2);
    yaw = rpy(3);

    halfYaw = yaw *.5;  
	halfPitch = pitch * .5;  
	halfRoll = roll *.5;
	
	cosYaw = cos(halfYaw);
	sinYaw = sin(halfYaw);
	cosPitch = cos(halfPitch);
	sinPitch = sin(halfPitch);
	cosRoll = cos(halfRoll);
	sinRoll = sin(halfRoll);
	
	 quat = [ cosRoll * cosPitch * cosYaw + sinRoll * sinPitch * sinYaw;
              sinRoll * cosPitch * cosYaw - cosRoll * sinPitch * sinYaw; 
			  cosRoll * sinPitch * cosYaw + sinRoll * cosPitch * sinYaw; 
			  cosRoll * cosPitch * sinYaw - sinRoll * sinPitch * cosYaw];

end

