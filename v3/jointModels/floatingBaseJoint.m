classdef floatingBaseJoint

    properties
        nv = {6}
        nq = {7}
        Xtree
    end
    
    methods        
        function [Xup, S, Sd, v] = jointKinematics(obj, q, qdot, vp)
            [XJ , S ] = jcalc( obj.jointType, q);   
            alc( obj.jointType, q);            
            Xup = XJ*obj.Xtree;
            if nargout > 2
                v  = Xup*vp + S*qdot;
                Sd = crm(v)*S;
            end
        end
    end
end

