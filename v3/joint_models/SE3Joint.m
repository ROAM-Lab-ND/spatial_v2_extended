classdef SE3Joint

    properties
        nv = {6}
        nq = {16}
        Xtree
    end
    
    methods        
        function [Xup, S, Sd, v] = kinematics(obj, q, qdot, vp)
            [XJ , S ] = jcalc( 'SE3', q);            
            Xup = XJ*obj.Xtree;
            if nargout > 2
                v  = Xup*vp + S*qdot;
                Sd = crm(v)*S;
            end
        end
    end
end

