classdef floatingBaseJoint

    properties
        nv = 6
        nq = 7
        bodies = 1
        output_body=1
    end
    
    methods        
        function [Xup, S, Sd, v] = kinematics(obj, Xtree, q, qdot, vp)
            [XJ , S ] = jcalc( 'Fb', q);              
            Xup = XJ*Xtree;
            if nargout > 2
                v  = Xup*vp + S*qdot;
                Sd = crm(v)*S;
            end
        end
    end
end

