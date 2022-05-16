classdef revoluteJoint

    properties
        nv = {1}
        nq = {1}
        axis
        Xtree
    end
    
    methods        
        function [Xup, S, Sd, v] = kinematics(obj, q, qdot, vp)
            [XJ , S ] = jcalc( obj.axis, q);            
            Xup = XJ*obj.Xtree;
            if nargout > 2
                v  = Xup*vp + S*qdot;
                Sd = crm(v)*S;
            end
        end
    end
end

