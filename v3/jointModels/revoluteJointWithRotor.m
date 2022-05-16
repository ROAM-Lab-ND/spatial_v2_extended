classdef revoluteJointWithRotor
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nv = {1}
        nq = {1}
        linkJointType
        rotorJointType
        Xtree
        gearRatio
    end
    
    methods        
        function [Xup, S, Sd, v] = jointKinematics(obj, q, qdot, vp)
            [XJ_link , S_link ] = jcalc( obj.linkJointType, q);
            [XJ_rotor, S_rotor ] = jcalc( obj.rotorJointType, q*obj.gearRatio);
            
            Xup = obj.Xtree;
            Xup(1:6 ,:) = XJ_rotor*Xup(1:6,:) ; 
            Xup(7:12,:) = XJ_link *Xup(7:12,:);

            S  = [S_rotor ; S_link];
            if nargout > 2
                v  = Xup*vp + S*qdot;
                Sd = crm(v)*S;
            end
        end
    end
end

