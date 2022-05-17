classdef revoluteJoint
    properties
        nv = 1
        nq = 1
        axis
        bodies = 1
        output_body=1
    end
    
    methods       
        function [Xup, S, Sd, v] = kinematics(obj, Xtree, q, qdot, vp)
            [XJ , S ] = jcalc( ['R' obj.axis], q);            
            Xup = XJ*Xtree;
            if nargout > 2
                v  = Xup*vp + S*qdot;
                Sd = crm(v)*S;
            end
        end        
    end
end

