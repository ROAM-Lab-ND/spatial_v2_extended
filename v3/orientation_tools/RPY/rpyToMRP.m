function [ mrp ] = rpyToMRP( rpy )
    mrp = rotToMRP(rpyToRot(rpy));
end

