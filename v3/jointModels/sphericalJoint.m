classdef sphericalJoint

    properties
        nv = {3}
        nq = {4}
        Xtree
    end
    
    methods        
        function [Xup, S, Sd, v] = jointKinematics(obj, q, qdot, vp)
            [XJ , S ] = jcalc( obj.jointType, q);            
            Xup = XJ*obj.Xtree;
            if nargout > 2
                v  = Xup*vp + S*qdot;
                Sd = crm(v)*S;
            end
        end
    end
end

