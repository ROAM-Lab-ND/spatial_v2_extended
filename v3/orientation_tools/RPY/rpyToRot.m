function [ R ] = rpyToRot( rpy )
% rpyToR takes earth-fixed sequenial roll-pitch-yaw euler angles to a
% rotation matrix
%   [R] = rpyToRot(rpy)
%   note: earth-fixed roll-pitch-yaw is the same as body-fixed
%   yaw-pitch-roll sequence 

    R = Rz(rpy(3))*Ry(rpy(2))*Rx(rpy(1));
end

