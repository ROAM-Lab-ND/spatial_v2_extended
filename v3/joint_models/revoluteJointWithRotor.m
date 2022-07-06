classdef revoluteJointWithRotor
    %revoluteJointWithRotor 
    
    properties
        nv = 1
        nq = 1
        jointAxis
        rotorAxis
        gearRatio
        output_body = 1
        bodies = 2
    end
    
    methods        
        function [Xup, S, Sd, v] = kinematics(obj, Xtree, q, qdot, vp)
            [XJ_link , S_link ] = jcalc( ['R' obj.jointAxis], q);
            [XJ_rotor, S_rotor ] = jcalc( ['R' obj.rotorAxis], q*obj.gearRatio);
            
            Xup = Xtree;
            
            Xup(1:6,:)   = XJ_link *Xup(1:6,:)  ;
            Xup(7:12 ,:) = XJ_rotor*Xup(7:12,:) ; 
            
            S  = [S_link ; S_rotor*obj.gearRatio];
            if nargout > 2
                v  = Xup*vp + S*qdot;
                Sd = crm(v)*S;
            end 
        end
        
    end
end

